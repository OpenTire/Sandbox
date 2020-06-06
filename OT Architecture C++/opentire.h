namespace opentire {

// Primary language-specific data types used throughout the OpenTire interface.

// Definitions below are pseudo-C++ but each type is intended to be a language-specific
// native type.  The type names themselves are also language-specific.  For example,
// in Matlab, a "vec3" is a 3-element column vector, e.g. zeros(3,1).  In Python a
// "vec3" is a 3-element column numpy vector.
//
// Distinct types below such as scalar, integer and bool, may be implemented
// with the same type in some languages (e.g. Matlab double).

// Each language will presumably provide some basic vector and matrix math functionality
// for vec2, vec3, and mat3: addition, subtraction, elementwise multiplication and division,
// matrix times vector, matrix or vector times scalar, matrix multiplication, etc.
//

// Class semantics notes:
// Types declared with "class" follow reference or "handle" semantics.
// two variables assiged the same handle object refer to the same object.
// The copy() method is supported for all handle objects.
// Any changes to a handle object will be visible to all references to that object.
//
// Types declared with struct (or using, belo) follow value semantics: assignments
// make a copy of the entire structure (so no copy() method is needed).
//
// Note that some classes are defined as classes within another class:
// For example, the tire_model input structure is declared nested within the
// tire_model class itself.  For languages that don't support this concept (e.g. Matlab),
// nested class names will be prepended by the class name of their containing
// class: the tire_model input class becomes tire_model_input.
//
// All classes and structures are intended to be extensible with the addition of new properties
// or methods, which will not affect existing code using the extended classes.
//
// Tire model implementations do not have to use every field in every structure;
// for example, the effect of tire inflation pressure input may not be implemented by all models.
// Models need not implement all of the tire model outputs, though this may restrict use of the
// model in some situations (e.g. dynamic simulations).
//

using scalar = double;    // 64-bit double-precision
using integer = int32_t;  // 32-bit signed integer
using bool = int32_t;     // an integer boolean
using string = std::string;

// Fixed-sized vectors.  Language may implement as
// small vectors.  Language provides appropriate operations such
// as multiplication, addition, dot and cross products, etc.
struct vec3 {
  scalar data[3];
};
struct vec2 {
  scalar data[2];
};
struct mat3 {
  scalar data[3][3];
};

// vector is a variable-sized, one-dimensional vector of scalars.
class vector {
  int size;
  scalar data[size];
};

// vectors of integers and strings of variable length
class integer_vector;
class string_vector;

// Type used for data that can take on different types
struct variant {};
class variant_vector;

// Two-dimensional matrix of scalars, with variable dimensions.
// Language provides implementation.
class matrix {};

//** OpenTire Classes

//** road_model base class for all road contact model implementations
class road_model {

  // compute_contact() input structure
  struct input {
    vec3 wheel_axis;
    vec3 wheel_center;
    scalar time;
    scalar distance;
  };

  struct contact {
    vec3 tire_contact;
    vec3 tire_contact_vel;
    // T.col(2) => wheel_axis
    mat3 T; // Tire transform: ground from local wheel/tire transform
    // W.col(3) -> road normal
    mat3 W; // Tydex-W transform: ground from W
    scalar road_grip_scale;
  };

  // Compute road contact information
  contact compute_contact(input input);
};

// Standard road model implementation for simple planar smooth road.
class plane_road_model : road_model {
  vec3 road_origin;
  vec3 road_normal;
  scalar road_grip_scale;
};

//*** tire_model : base class for all tire model implementations.

// Coordinate System and Units.
//
// All input and output values throughout OpenTire are in SI units (mm, N, kg, rad)
// and use the ISO/Tydex-W coordinate system. ISO is a right-handed
// system with X+ forward, Y+ to the left, and Z+ up. Rotations follow
// the right-hand rule: grasp an axis with the right hand with thumb
// pointing in positive direction: the curl of the hand defines the
// positive rotation direction.
// For convenience a function is provided for transforming quantities
// from ISO to and from other coordinate systems.

class tire_model {

  struct options {
    int log_level = 0; // Logging level, 0 = none.
    // Non-zero values produce model-dependent output in output.log_data.

    // Any other common options?  We can define common options that aren't used by all models.
    // symmetric?
    // use_turn_slip?
  };

  // Tire model scale factors.
  // These are multiplied by scale factors of tire model itself (if any),
  // and road_model output, and are specified as input to
  // compute_steady() and compute_dynamic() functions.
  // Not all models support all scale factors.

  struct input_scale_factors {
    // Intent is to support the 7 scale factors supported by TNO.
    // At time of writing I'm not sure what they are: can't find the docs!
    vec2 grip_adjust; // lon and lat grip adjustment
    scalar velocity_grip_adjust; // speed-sensitive grip adjustment
    // others TBD?
  };

  // Common, extensible tire output structure
  // This structure is never created directly by client programs;
  // instead it is returned by one of the compute_*() functions.
  // This allows models to extend this structure in the future.
  struct output {
    scalar right_tire = 0;

    vec3 force; // [N] Fx, Fy, Fz
    vec3 moment; // [N*m] Mx, My, Mz

    vec3 tire_contact; // [m] Tire contact point
    vec3 tire_contact_vel; // [m/s] Velocities of tire contact point
    vec2 tire_slip_vel; // [m/s] Tire longitudinal and lateral slip velocities

    scalar wheel_omega; // [rad/s] wheel angular velocity
    scalar rolling_radius; // [m] effective rolling radius

    scalar kappa; // [-1..1] Longitudinal slip
    scalar alpha; // [rad] Lateral slip
    scalar gamma; // [rad] Camber/inclination angle
    scalar path_curvature; // [1/m] Reciprocal of path radius 1/R
    scalar phi; // [1/m] effective turnslip (path curvature plus gamma-induced)

    scalar loaded_radius; // [m]
    scalar rho; // [m]
    scalar vertical_stiffness; // [N/m]
    scalar vertical_damping; // [N/(m/s)]

    scalar trail; // [m] longitudinal and lateral pneumatic trails
    scalar grip; // [-] lon and lat grip

    scalar relaxation; // [m] lon and lat relaxation lengths
    scalar slip_stiffness; // nominal slip stiffnesses (dFx/dkappa, dFy/dalpha, dFy/dgamma, dMz/dphi)
    scalar contact_size; // [m] X and Y contact patch dimensions (full size, not half)

    mat3 W; // Rotation to inertial frame from Tydex W frame
    mat3 R; // Rotation to inertial frame from wheel frame

    vector states_dot; // State time derivatives

    vector algebraic_loops; // Internal algebraic loop variables

    vector log_data; // Additional output data, controlled by options.log_level
  };

  // Base class for compute_steady() input
  // This structure is never created directly by client programs;
  // instead it is created by a tire object and returned.  This allows
  // models to extend this structure in the future.
  struct steady_input {
    scalar right_tire; // [bool] 1 if tire mounted on right, 0 left.

    scalar Fz = nan; // [N] Vertical force in N (0 if not known)
    scalar loaded_radius = nan; // [m] loaded_radius (nan/0 if not known)
    scalar path_curvature = nan; // [1/m] 1/R
    scalar Vx = nan; // [m/s]
    scalar wheel_omega = nan; // [rad/s] (nan/0 if not known)
    scalar kappa = nan; // [-] (nan/0 if not known)
    scalar alpha = nan; // [rad]
    scalar gamma = nan; // [rad]
    scalar pressure = 0; // [pascal]  (0 for nominal/unknown)

    input_scale_factors scale_factors;
  };

  struct dynamic_input {
    scalar right_tire;    // [bool] 1 if mounted on right, 0 if left

    integer solver_state; // [-] Code indicating state of solver.
                          // Similar to STI JOBFLG -- to be defined.

    scalar time;          // [s] time since simulation start (optional)
    scalar distance;      // [m] simulated distance traveled (optional)

    vec3 frame_omega;     // [rad/s] angular velocity of rotating frame

    vec3 wheel_center;    // [m] wheel center in inertial frame
    vec3 wheel_center_vel; // [m/s] wheel center velocity in inertial frame

    vec3 wheel_axis;      // [-] wheel axis normal in inertial frame

    scalar wheel_angle;   // [rad] wheel rotation angle about wheel_axis (nan if not known)
    scalar wheel_omega;   // [rad/s] wheel angular velocity about wheel_axis

    vec2 grip_adjust;     // [-] grip adjustment
    scalar pressure;      // [pascal] inflation pressure

    input_scale_factors scale_factors;

    vector states;  // Differential equation states: these states are integrated by the calling simulation environment,
                    // using the time derivatives states_dot returned by the model.

    vector algebraic_loops; // internal algebraic loop variables: output variables that are also inputs.
                            // algebraic_loops (if any) are model-specific.
  };

  //*** Public tire_model properties
  // (These may be implemented with property get/set functions,
  // but they are accessed as normal data properties)
  // Default implementation provided.

  //*** Public tire_model methods

  // Get model version number, formatted as:
  // major * 100 + minor * 10 + revision
  integer model_version;

  // Filename (may be empty)
  string filename;

  // Comment string (multiline format, e.g. ADAMS-specific header info in MF 6.2 TIR files)
  string comment;

  // true if tire object is valid, false if load() failed, etc.
  bool is_valid;

  // Return a new tire_model_options with default settings.
  options new_options();

  // Return a new steady_input structure with default values
  steady_input new_steady_input(bool right_tire, road_model road);

  // Return a new dynamic_input structure with default values
  dynamic_input new_dynamic_input(bool right_tire, road_model road);

  // Load parameters from a file
  bool load(string filename);

  // format 0 - standard text file format
  // format 1 - standard binary format (classname followed by parameter vectors in binary)
  // formats 2-99 are reserved.
  // Other model-dependent formats start at 100.
  bool save(string filename, int format);

  // Get and set parameters by index
  vector get_parameter_by_index(int index);
  bool set_parameter_by_index(int index, vector parameters);

  int get_parameter_index_by_name(string names);

  // Common parameter acces helper functions
  vector get_parameter_vector();
  bool set_parameter_vector(vector params);

  // Vector variations
  vector get_parameter_by_index(int_vector indices);
  bool set_parameter_by_index(int_vector indices, vector params);

  vector get_parameters_by_name(string_vector names);
  bool set_parameter_by_index(string_vector names, vector params);

  //*** Compatibility functions

  // TNO dfeval() (and MFeval) steady-state computation function.
  // input is an N-by-6 matrix of input values;
  // output is an N-by-40 matrix of output values.
  // See TNO dfeval.m for detailed description and
  // definition of the matrix columns.
  matrix dfeval(matrix input, options options);

  //*** Analysis functions
  // These functions may be best implemented as functions taking
  // tire_model objects as arguments, rather than as methods as described here.
  // These operations will tend to be very platform-specific: producing charts
  // and writing output to a console window.  By providing them as functions
  // they can be easily adapted to the underlying language and platform.

  // Prints a comparison of the properties of this tire and other_tire.
  bool diff(tire_model other_tire);

  // Arguments TBD.  Will support specifying up to 3 named inputs and their ranges,
  // and one or two output axes.
  bool plot(x_inputs, y_inputs, z_inputs, output_x, output_axis_1, output_axis_2);

  //*** Class Implementation Methods
protected: // only available to tire model implementors
  // A model implementation is free to override compute_steady() and compute_dynamic()
  // and implement them directly, but it is often convenient to have a single implementation
  // function for both.
  output compute_common(steady_input, dynamic_input, road_model road, options options);

  // By providing an implementation of get_parameter_info(), open_tire will provide
  // a default implementation of many of the methods above, such as: load, save
  // diff, get/set_parameter_by_name/index, and more.
  //
  struct parameter_info {
    string name;
    string description;
    integer type;
    integer parameter_index;
    scalar supported_version;
    scalar default_value;
  };
  // vector of parameter_infos
  struct parameter_info_vector {};

  // return table of info used to implement various standard methods
  parameter_info_vector get_parameter_info();
};  // end of class tire_model

//** Tire model implementation subclasses
// (of course, final version of opentire may or may not include all of these models)

// simple_tire - A minimal tire model consisting by a minimum set of parameters, intended to be
// used as the basis for new models.  Can also be used in its own right for simple studies.
class simple_tire : tire_model;

//** mmftire - Magic Formula implementation as described in Chapter 4 of Tire and Vehicle Dynamics, by H.B.Pacejka.
// and implemented by TNO.  Supports versions 5.2, 6.0, 6.1, 6.2; will support future versions as they are defined.
// A goal is to implement the full MF-SWIFT model in the future.
class mftire : tire_model;

//** brush_model - Brush model as described in Chapter 3 of Tire and Vehicle Dynamics, by H.B.Pacejka.
class brush_model : tire_model;

//** TYR501 - Empirical model described in Multibody Systems Approach to Vehicle Dynamics by D. Harty.
// Another good starting point for new model development.
class TYR501 : tire_model;

//** Standard functions

// Change an object (vector, matrix, tire_input, etc)
// to/from the standard ISO/Tydex-W frame to another coordinate system.
// This function is its own inverse: the same call is used to convert from ISO to another system
// as well as the other system back to ISO.
//
// This function can be used by clients as well as within model implementatoins.
//
// "object" parameter can be any of:
//      vec3
//      mat3
//      tire_model::dynamic_input
//      tire_model::steady_input
//      tire_model::output
//

// Standard names for supported coordinate system identifiers
// (Adapted_SAE is used in Tire and Vehicle Dynamics by Pacejka)
enum coordinate_system_id : integer {
  ISO = 0,
  SAE = 1,
  Adapted_SAE = 2
};

variant change_coordinate_system(variant object, coordinate_system_id other_system);

} // namespace opentire
