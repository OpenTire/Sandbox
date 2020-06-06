%% OpenTire class interface proposal
%
% Abstract tire model class

% All tire models derive from opentire.tire(), and have reference (handle) semantics.
classdef tire < handle
    
    function ok = load(filename)
    end
    
    function ok = save(filename)
    end
    
    function options = new_options()
    end
    
    % Steady state methods
    function input = new_steady_input(right_tire)
    end
    
    %% compute_steady method
    % Computes steady-state forces, moments, and other tire_output values.
    % Some tire_input fields are optional, indicated by the value nan.
    % Computed by compute_steady() if not present (e.g. typically loaded_radius, wheel_omega)
    % Use new_* functions to obtain instances of input and options structs.
    function output = tire.compute_steady(input, options)
    end
    
    % Dynamic methods
    function [input,states,algebraic_loops] = new_dynamic_input(right_tire, options)
    end
    
    %% compute_dynamic method
    % algebraic_loops is [] or a vector of algebraic-loop values (outputs that are also inputs).
    % Integrators should call this function with values from previous time step.
    % Equilibrium solvers should solve for (output algebraic_loops - input_algebraic_loops) == 0.
    % states is [] or a vector of states whose time-derivatives are computed by the tire model.
    % state time derivatives are returned in state_dots.
    
    function [output,state_dots,algebraic_loops] = compute_dynamic(time, input, road, states, algebraic_loops)
    end
    
end % classdef tire

% The two tire classes provided in initial version of OpenTire

% Magic Formula 6.2 tire model
classdef mftire < opentire.tire
    % implementation of all tire functions
end

% Brush model as described in TVD Chapter 3.
classdef brush_tire < opentire.tire
    % implementation of all tire functions
end

% Create/copy/load/save examples
% tire = opentire.mftire() % create an empty tire (can't yet be used for computation)
% tire = opentire.mftire('folder/filename.tir'); % create a tire from the file folder/filename.tir
% right_tire = copy(left_tire); % create a copy of a tire.
% tire.load('folder/filename.tir');
% tire.save('folder/secondname.tir');
% tire.save();


%% Model computation options structure.
% Each class has its own options subclass, which cannot be created
% directly: instead use tire.get_default_options() to create an instance.
%
classdef tire_options
    % no common options I can think of
    % but if there were, they'd go here
end

function options = tire.new_options()
end

%% Steady state computation function
function output = tire.compute_steady(input, options)
end

% Steady-state input base class
% Tire model classes can extend this by subclassing.
% Tire models that don't support some of the following parameters
% can simply ignore them (e.g. pressure, curvature)
%
% NOTE: This structure is not created directly:
% create one with the tire.get_steady_input() function.
%
classdef steady_input
    right_tire  % [bool] 1 if tire mounted on right, 0 left.
    Fz          % [N] Vertical force in N (0 if not known)
    loaded_radius % [m] loaded_radius (nan/0 if not known)
    curvature   % [1/m] path curvature, 1/R
    Vc          % [m/s]
    wheel_omega % [rad/s] (nan/0 if not known)
    kappa       % [-] (nan/0 if not known)
    alpha       % [rad]
    gamma       % [rad]
    pressure    % [pascal]  (0 for nominal/unknown)
    % Tire models can extend this class via subclassing, to add, e.g.:
    % temperature(3) % [C] inner, middle, outer tire temps
    % wear(3)        % [-] tire wear inner, middle outer
end

%% new_input method
function input = tire.new_input(right_tire)
end

% mftire class options example
classdef mftire_options < tire_options
    turn_slip_model  % enable turn_slip model
    transient_model  % selects transient model to use 0 = none, 1 = restricted nonlinear, etc
    no_wheel_lift    % Allow negative FZ if wheel leaves the ground: useful for equilibrium solving
    use_alpha_gamma_star % use sin(alpha) instead of alpha throughout MF computations
end

% Both compute_steady_state() and compute_dynamic() produce the following output structure
classdef tire_output
    log_data        % additional logging output data enabled by options (or [])
    
    force(3)    % Fx, Fy, Fz  as a 3-vector
    moment(3)   % Mx, My, Mz  as a 3-vector
    
    wheel_omega
    rolling_radius
    
    phi         % [1/m] Effective turnslip (path curvature plus gamma-induced)
    
    loaded_radius   % [m]
    vertical_stiffness [N/m]
    vertical_damping % [N/(m/s)]
    
    trail(2)        % [m] longitudinal and lateral pneumatic trails
    grip(2)         % [-] lon and lat grip
    contact_size(2) % [m] X and Y contact patch dimensions (full size, not half)
    relaxation(2)   % [m] lon and lat relaxation lengths
    W(3,3)          % Rotation from Tydex W frame to vehicle frame
    R(3,3)          % Rotation from wheel frame to vehicle frame
end

%% change_coordinate_system method
% All input and output classes have the change_coordinate_system
% method used to transform input and output structures
% from one coordinate system to another.
% Note that compute_steady() and compute_dynamic() must
% be called with inputs in the ISO coordinate system,
% and will always return outputs in ISO as well.  It is
% the responsibility of the caller to use
% change_coordinate_system when needed.
ISO_SYSTEM = 0;
SAE_SYSTEM = 1;
ADAPTED_SAE_SYSTEM = 2;

function transformed = change_coordinate_system(input_or_output, from, to)
end

%% abstract input class for compute_dynamic()
classdef dynamic_input
    right_tire        % [bool] 1 if mounted on right, 0 if left
    
    % Angular velocity of rotating frame
    frame_omega(3)
    wheel_center(3)   % [m] wheel center in inertial frame
    wheel_center_vel(3) % [m/s] wheel center velocity in inertial frame
    wheel_axis(3)     % wheel axis normal in inertial frame
    wheel_omega       % [rad/s] wheel angular velocity about wheel_axis
    
    grip_adjust(2)    % grip adjustment
    pressure          % inflation pressure
end

%% Abstract road_model interface.
classdef road_model
    
    %% Road model single-contact-point function.
    % returns contact point, surface normal,
    % and surface curvature (second partial derivatives)
    % all in the inertial frame.
    % distance parameter is integral of vehicle Vx and may be 0.
    % Modeler must arrange to integrate vehicle Vx and pass this
    % parameter to road models that support it.
    function [point,normal,grip_adjust,dzzdx,dzzdy] = road.compute_contact(wheel_axis, wheel_center, distance)
    end
    
end

%% Simple single-contact-point model implementation
%
% Simple planar road model
classdef plane_road_model < opentire.road_model
    origin(3)   % plane origin in inertial frame
    normal(3)   % plane normal in inertial frame
    grip_adjust % grip adjustment
end
