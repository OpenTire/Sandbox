%% opentire.mftire version 1.1.0 - opentire.mftire class
%
% This class implements the Magic Formula tire model as described in
% "Tire and Vehicle Dynamics, 3rd Edition", by Hans Pacejka.
%
% Based on MFeval, written by Marco Furlan

classdef (Sealed = true) mftire < opentire.tire

%%
%
% References:
% [1] Title:    Tire and Vehicle Dynamics
% Author:       Hans Pacejka
% Edition:      3, revised
% Publisher:    Elsevier, 2012
% ISBN:         0080970176, 9780080970172
% Length:       672 pages
%
% [2] the *latest revision* of the following paper:
%   Besselink, I. J. M. , Schmeitz, A. J. C. and Pacejka, H. B.(2010) 'An
%   improved Magic Formula/Swift tyre model that can handle inflation
%   pressure changes', Vehicle System Dynamics
%
% The published version of this paper is slightly out of date;
% the latest revision can be found at:
%
%   https://pure.tue.nl/ws/files/3139488/677330157969510.pdf
%
% [3] MSC.Software "Using the PAC2002Tire Model"
% which can be found at:
%
%   http://mech.unibg.it/~lorenzi/VD%26S/Matlab/Tire/tire_models_pac2002.pdf
%
% The loaded radius model is implemented from:
%
% [4] R.H.M.T. van der Hofstad (2010)  Master's Thesis
% 'Study on improving the MFSwift tyre model'
%
% which can be found at:
%
%   http://www.mate.tue.nl/mate/pdfs/11495.pdf
%
% This code is based on MFeval by Marco Furlan
% Script version:        1.5, 4.0
% Code developed by:     Marco Furlan
% Email address 1:       mfurlant@jaguarlandrover.com
% Email address 2:       marcofurlan92@gmail.com

% Coding conventions to assist with translation to C++
% 1) If statement conditionals are always parenthesized
% 2) If statements have '%then' comment after parenthesized conditional
% 3) Constants such as 1, 0, -1 are passed to functions as 1.0, 0.0, etc
% 4) Fractions such as 1/3 are written as 1.0/3.0
% 5) POW() function is used instead of ^ operator
% 6) SQR() function is used instead of ^2
% These conventions make translating this code to C++ much easier.

methods(Access = public)

%% Constructor
% Construct from tire filename or param vector
function tire = mftire(params)
  if nargin == 0
    tire.is_valid = true;
    tire.on_params_changed();
    
  elseif (ischar(params)) %then
    
    % string: load from filename
    tire.load(params);
    
  elseif (isvector(params) && length(params) > 1) %then
    
    % initialize from a parameter vector
    tire.set_param_vector(params);
    
  elseif (params == 0) %then
    % return a copy with reference parameters
    tire.is_valid = true;
    tire.on_params_changed();
    
  else
    warning('opentire.mftire: Unsupported constructor parameter');
    tire = []; % return empty handle to indicate error
  end
end

%% Create options structure
function options = new_options(tire)
  options.want_turn_slip = false;
  options.transient_model = 0;
  options.use_alpha_gamma_star = false;
  options.use_adapted_SAE = false;
  options.no_wheel_lift = false;
end

%% Create input structure with nominal values for steady-state simulations
function input = new_steady_input(tire, right_tire)
  
  if (nargin < 2)  %then
    right_tire = false;
  end
  
  input = opentire.steady_input();
  input.right_tire = right_tire;
  input.Fz = tire.FNOMIN;
  input.Vx = tire.LONGVL;
  input.loaded_radius = nan;
  input.kappa = 0;
  input.alpha = 0;
  input.gamma = 0;
  input.path_curvature = 0;
  input.pressure = tire.INFLPRES;
end

function input = new_dynamic_input(tire, right_tire)
  if (nargin < 2)  %then
    right_tire = false;
  end
  
  input = opentire.dynamic_input();
  input.right_tire = right_tire;
end


%% compute_steady, compute_dynamic -- Tire model computation functions
%
% These functions computes tire forces and moments from the supplied input.
%
% Two kinds of computations are supported: steady state and dynamic.
%
% Steady state computations simulate a tire test bench, where
% the key tire parameters such as vertical load, speed, camber, 
% slip angle and longitudinal slip are fixed, and the resulting
% forces and moment are computed.
%
% Dynamic computations are for incorporation in a multibody dynamic
% simulation.  The state of the tire is specified as positions
% and rotations in an inertial space along with their time derivatives,
% from which the resulting forces and moments are computed.
%
% The input structure has separate steady-state and dynamic input fields.
% Steady state computations are performed when either the input.Fz
% or input.loaded_radius fields are not Nan -- in that case,
% the other state variables kappa, alpha, gamma must be specified.
% If loaded_radius is NaN, then the simulation is done at the specified
% Fz vertical load.  Otherwise, loaded_radius is considered fixed,
% and the resulting vertical load at loaded_radius is computed.
% (Because fixed loaded_radius simulations are done iteratively,
% it is helpful to provide an input.Fz value as an initial estimate
% of the load at the given loaded_radius).
%
% Dynamic simulations are performed when both input.Fz and input.loaded_radius
% are NaN.  The dynamic input members must all be specified; the
% tire slips, camber, and loaded radius and vertical load are computed
% along with the tire forces and moments.
%
% Both steady-state and dynamic computations share the
% road and road_grip inputs.  The road plane is defined by
% by an origin and a normal to the road surface.  The road
% grip inputs allow for changing lateral and longitudinal
% grip as a simulation input (a function of time, or position, etc).
%
% Algebraic loops
%
% As it turns out, the correct computation of forces and moments
% may depend on the very forces being computed -- in other words,
% the outputs can also be inputs.
%
% The vertical stiffness of the tire can sometimes be affected
% be affected by the horizontal and lateral forces on the tire.
% The stiffness affects the vertical load produced by a given
% loaded radius, which in turn affects the lateral and horizontal
% tire forces. This coupling of (Fx,Fy)->stiffness->Fz->(Fx,Fy)
% is known as an algebraic loop.
%
% In the presence of an algebraic loop, the computed forces
% aren't correct until the input forces are equal to the output forces.
% Force is included as a dynamic input for this reason (although
% only force(1) and force(2) are used).
%
% Not all tires have an fx,fy algebraic loop: the 
% has_algebraic_loop() function indicates whether a tire
% has an algebraic loop, and also which forces are coupled.
%
% In practice, for a dynamic simulation-in-time, it is usually
% sufficient to set the input forces to the output forces from
% a previous time step. But sometimes no reasonably accurate
% estimate is available (such as during equilibrium solving)
% and the algebraic loop must be closed by iteratively solving
% for output forces == input forces.  The function
% close_algebraic_loop() will do this for you, at a cost of
% 3 or 4 extra evaluations.
%
% TVD book conventions are "adapted SAE", but model uses ISO conventions: See Appendix 1.
% Differences are:
%   * Not a right handed system
%   * Sign of Y axis (Fy, Vc(1), Vs(1))
%   * Sign of Z axis rotations (Mz, psi/psi_dot)
%   * alpha and phi/phi_t are negated, but kappa has same sign.
%   * positive camber rotations counter clockwise (viewed from rear) in both systems, but:
%   * positive camber is towards inside of curve rather than outside as in ISO
%   * omega_wheel = -omega rather than omega_wheel == omega as in ISO
%  
% Differences are sign of Y axis (Fy, Vc(1), Vs(1)) and in Z axis rotations (Mz, psi_dot)
% Turn slip signs are the same: sign(phi) = -sign(psi_dot) in adapted SAE,
% sign(phi) = sign(psi_dot) in ISO.

function [output, states_dot, algebraic_loops] = compute_dynamic(t, input, road, states, algebraic_loops, options)
  [output, states_dot, algebraic_loops] = t.compute_common([], input, road, states, algebraic_loops, options);
end

function output = compute_steady(t, input, options)
  [output,~,~] = t.compute_common(input, [], [], [], [], options);
end

function [output, states_dot, algebraic_loops] = compute_common(t, steady_input, dynamic_input, road, states, algebraic_loops, options)
%% Process inputs
  output = opentire.tire_output;
  
  states_dot = nan(size(states));

  if (~t.is_valid) %then
    fprintf('opentire.mftire %s: Invalid tire params\n', t.filename);
    return;
  end
  
  if (isempty(dynamic_input)) %then
    steady_state = true;
  
    if (isfinite(steady_input.loaded_radius))
      % If input.loaded_radius is not NaN, then
      % we want to compute Fz at a fixed loaded radius.
      % This is done iteratively (i.e. slowly).
      % In this case, input.Fz is a hint used to speed up the iteration;
      % if not specified (i.e. NaN) then t.FZNOMIN is used as an
      % iteration starting point.
      %
      if (isfinite(steady_input.loaded_radius)) %then
        % iteratively solve for Fz at given loaded radius
        output = t.fixed_loaded_radius_eval(steady_input, output, options);
        return;
      end
    end
    
    right_tire  = (steady_input.right_tire ~= 0);
    pressure = steady_input.pressure;
    
  else
    steady_state = false;
    
    right_tire  = (dynamic_input.right_tire ~= 0);
    pressure = dynamic_input.pressure;
    
  end

  want_turn_slip  = options.want_turn_slip;
  transient_model = options.transient_model;
  use_alpha_gamma_star = options.use_alpha_gamma_star;
  use_adapted_SAE = options.use_adapted_SAE;
  no_wheel_lift   = options.no_wheel_lift;
  
  if (~isfinite(pressure) || pressure == 0) %then
    pressure = t.INFLPRES;
  end
  pressure = CLAMP(pressure, t.PRESMIN, t.PRESMAX);
  
  % Normalized inflation pressure
  pi0 = t.NOMPRES; % Reference pressure
  dpi = 0;
  if (pi0 > 0) %then
    dpi = (pressure - pi0)/pi0; % [Eqn (4.E2b) Page 177 - Book]
  end
  
  % mirror data input and output if on opposite side
  % Map input state into values consistent with the mounted tire side
  side_flip = 1;
  if (right_tire ~= (t.TYRESIDE ~= 0)) %then
    side_flip = -1;
  end
  
  input_Fz = nan;
  tire_contact_vel = [nan;nan;nan];
  rolling_radius = nan;
  
  R0 = t.UNLOADED_RADIUS;
  V0 = t.LONGVL; % Nominal speed (>0)
  inv_V0 = 1/V0;
  
  % Nominal stiffness (pressure corrected)
  CFz0     = t.VERTICAL_STIFFNESS;
  CFz = CFz0 * (1 + t.PFZ1*dpi);  % (Paper Eq. 5) - Vertical stiffness adapted for tire inflation pressure

  if (steady_state) %then
    
    % IMPORTANT NOTE: Vx = Vcx [Eqn (7.4) Page 331 - Book]
    input_Fz = steady_input.Fz;
    input_Vx = steady_input.Vx;
    
    input_kappa = steady_input.kappa;
    input_alpha = steady_input.alpha;
    gamma = steady_input.gamma;

    % compute wheel_axis from alpha and gamma
    [sin_input_alpha, sin_gamma] = two_sin(input_alpha, gamma);
    [cos_input_alpha, cos_gamma] = two_cos(input_alpha, gamma);
    wheel_axis = [sin_input_alpha*cos_gamma; cos_input_alpha*cos_gamma; sin_gamma];

    gamma = gamma * side_flip;
    sin_gamma = sin_gamma * side_flip;

    wheel_omega = steady_input.wheel_omega;
    if (~isfinite(wheel_omega)) %then
      assert(isfinite(input_kappa));
      [rolling_radius, wheel_omega] = t.compute_rolling_radius_and_wheel_omega(input_Fz, input_kappa, input_Vx, CFz);
    end
    
    tire_contact_vel = [input_Vx/cos_input_alpha; 0; 0];

    % 1/path_radius about road surface normal which points upward.
    % Positive curvature values indicate rotation to the left.
    phi_t = steady_input.path_curvature;
    if (use_adapted_SAE) %then
      phi_t = -phi_t;
    end
    % Disable turn slip calculations if curvature is not defined
    if (~isfinite(phi_t)) %then
      want_turn_slip = false;
      phi_t = 0;
    end
    
    frame_omega = [0;0;phi_t * input_Vx];
    
    road_normal = [0;0;1];

    % R = transform from tire frame to vehicle frame
    % R(:,2) = wheel rotation axis normal

    R = zeros(3,3);
    R(:,1) = vec3_normalize(vec3_cross(wheel_axis, road_normal));
    R(:,2) = wheel_axis;
    R(:,3) = vec3_cross(R(:,1), wheel_axis);

    % W = transform from Tydex W tire frame to vehicle frame
    % W(:,3) == road surface normal, pointing upwards
    % 
    W = zeros(3,3);
    % Z axis: road surface normal, pointing upwards.
    W(:,3) = road_normal;
    % X axis: intersection of wheel plane with ground plane
    W(:,1) = vec3_normalize(vec3_cross(wheel_axis, W(:,3))); % R(:,2) = wheel_axis
    % Y axis: projection of wheel axis onto ground plane
    W(:,2) = vec3_cross(road_normal, W(:,1));
    
    road_grip = 1;
    
  else % compute_dynamic()

    wheel_axis = dynamic_input.wheel_axis;
    wheel_omega = dynamic_input.wheel_omega;
    frame_omega = dynamic_input.frame_omega;
    wheel_center = dynamic_input.wheel_center;
    wheel_center_vel = dynamic_input.wheel_center_vel;
    
    time = dynamic_input.time;
    distance = dynamic_input.distance;
    
    [tire_contact, road_curvature, R, W, road_grip] = road.compute_contact(wheel_axis, wheel_center, time, distance);

    radial = tire_contact - wheel_center;
    loaded_radius = vec3_norm(radial);

    % W' * R = rotate_x(gamma) = [1 0 0;
    %                             0 cos(gamma) -sin(gamma);
    %                             0 sin(gamma) cos(gamma)]:
    % so gamma can be extracted from elements (3,2:3) of W'*R.
    % (we use atan() instead of asin or acos so
    % we can use our fast estimated atan function)
    sin_gamma = vec3_dot(R(:,2), W(3,:)) * side_flip;
    cos_gamma = vec3_dot(R(:,3), W(3,:));
    gamma = atan(sin_gamma / cos_gamma);

    % Contact point velocity in W frame x and y axes
    % IMPORTANT NOTE: Vx = Vcx [Eqn (7.4) Page 331 - Book]
    % It is assumed that the difference between the wheel center longitudinal
    % velocity Vx and the longitudinal velocity Vcx of the contact center is
    % negligible
    tire_contact_vel = wheel_center_vel + vec3_cross(frame_omega, radial);
  end
  
  Vc = (tire_contact_vel' * W)'; % map to W frame
  Vc(2) = Vc(2) * side_flip;
  sign_Vcx = SGN(Vc(1));
  Vcx_prime = make_nonzero(Vc(1));

  % Compute phi_t as 1/path_radius in the W ground plane
  % Positive phi_t indicate rotation to the left, when viewed from above.
  if (~steady_state) %then
  end

  stiffness = nan;
  rho = nan;
  n_omega = wheel_omega * R0*inv_V0;
  R_omega = R0 * (t.Q_RE0 + t.Q_V1 * SQR(n_omega));
  
  if (~steady_state) %then
    
    omega_W_z = vec3_dot(frame_omega, W(:,3));
    phi_t = omega_W_z / Vcx_prime;
    
    input_Fx = 0;
    input_Fy = 0;
    if (~isempty(algebraic_loops))
      input_Fx = algebraic_loops(1);
      input_Fy = algebraic_loops(2) * side_flip;
    end
  
    %% Compute load from loaded radius
    % Compute vertical force from loaded radius and possibly previous lat & lon loads
    % Radius functions use unclamped slip, gamma, and pressure values.
    %function [input_Fz, stiffness, rho] = load_from_loaded_radius(t, loaded_radius, R_omega, n_omega, input_Fx, input_Fy, gamma, dpi, no_wheel_lift)
    
    % Handle rim bottoming
    % REVIEW: The wheel bottoming model described starting on Page 11 of [3]
    % is more elaborate than this. It also states that the additional Fz due to
    % bottoming should be included in Mx calcs, but not Fx,Fy,My and Mz.
    % So this code is disabled for now.
    bottom_stiffness = 0;
    bottom_Fz = 0;

    % Modeled as a parallel spring
    % REVIEW: but maybe more correct as a spring
    % that changes stiffness to BOTTOM_STIFF past min_loaded_radius.
    %   bottom_offset    = t.BOTTOM_OFFST;
    %   if (bottom_offset ~= 0) %then
    %     min_loaded_radius = t.RIM_RADIUS + bottom_offset;
    %     if (loaded_radius < min_loaded_radius) %then
    %       
    %       bottom_stiffness = t.BOTTOM_STIFF;
    %       % a reasonable default for bottom stiffness if not specified
    %       if (bottom_stiffness == 0) %then
    %         bottom_stiffness = t.VERTICAL_STIFFNESS * 10;
    %       end
    %       
    %       bottom_Fz = (min_loaded_radius - loaded_radius) * bottom_stiffness;
    %     end
    %   end

    [KR1, KR2, rho] = compute_rho_quadratic(t, R_omega, n_omega, input_Fx, input_Fy, gamma, dpi, no_wheel_lift, loaded_radius);

    input_Fz = KR1*rho + KR2*SQR(rho) + bottom_Fz;

    % treat as parallel springs
    stiffness = (KR1 + 2*KR2*rho) + bottom_stiffness;
    
    %end 
    %[input_Fz, stiffness, rho] = t.load_from_loaded_radius(loaded_radius, R_omega, n_omega, ...
    %    input_Fx, input_Fy, gamma, dpi, options.no_wheel_lift);
  end
  
  unclamped_Fz = input_Fz; % vertical load

  low_load_scale = 1;
  if (unclamped_Fz < t.FZMIN) %then
    low_load_scale = unclamped_Fz * (unclamped_Fz / t.FZMIN);
  end
  
  Fz = CLAMP(unclamped_Fz, t.FZMIN, t.FZMAX);
  
  % Normalized nominal load
  Fz0 = t.FNOMIN;
  inv_Fz0 = 1/Fz0;
  Fz0_prime =  t.LFZO*Fz0; % [Eqn (4.E1) Page 177 - Book]
  inv_Fz0_prime = 1/Fz0_prime;

  % Normalized vertical load
  nFz = Fz*inv_Fz0;
  unclamped_nFz = unclamped_Fz*inv_Fz0;
  
  nFz_prime = Fz*inv_Fz0_prime;
  dFz = nFz_prime - 1; % == (Fz - Fz0_prime)*inv_Fz0_prime; % [Eqn (4.E2a) Page 177 - Book]

  % Low speed scaling
  low_speed_scale = 1;
  if (abs(Vc(1)) < t.VXLOW) %then
    low_speed_scale = SMOOTHSTEP( abs(Vc(1)) / t.VXLOW );
  end
  
  epsilon_v = 1e-6;
  Vc_prime = sqrt(SQR(Vc(1)) + SQR(Vc(2))) + epsilon_v; % [Eqn (4.E6a) Page 178 - Book]
  cos_alpha_prime = Vc(1) / Vc_prime; % [Eqn (4.E6) Page 177 - Book]

  inv_Vcx_prime = 1 / Vcx_prime;
  
  gamma_star = sin_gamma;
  if (~options.use_alpha_gamma_star) %then
    gamma_star = gamma;
  end
  
  unclamped_gamma = gamma;
  gamma = CLAMP(unclamped_gamma, t.CAMMIN, t.CAMMAX);

  %% Spin slip parameter
  
  %function [phi, nphi, eps_gamma, zeta3] = compute_spin_slip(t, phi_t, sin_gamma, wheel_omega, inv_Vcx_prime, dFz, R0, want_turn_slip)

  phi = 0;
  nphi = 0;
  eps_gamma = 0;
  zeta3 = 1;
  
  phi_t = phi_t * sign_Vcx;
  phi_t = phi_t * side_flip;

  if (want_turn_slip) %then

    % Spin slip is defined in Eqs. 2.15-2.18 of the book, in the "adapted SAE" system.
    % Fig. 2.6 on page 66 shows that camber increases the magnitude of 
    % turnslip when the tire is cambered towards the center of the curve.
    % A tire in a right turn as shown in 2.6 will have negative path curvature
    % and a positive camber angle. By the right hand rule in the ISO system,
    % a positive wheel rotation with positive camber will induce a negative angular
    % velocity about the Z axis:
    % 
    %   psi_dot_omega = wheel_omega * -sin_gamma;
    % 
    % This corresponds to Eq. (7) of [3], which shows the negative sign
    % in front of the (1-eps_gamma) term.

    % IMPORTANT NOTE:
    % phi and phi_t above are computed using Vcx_prime,
    % rather than Vc_prime as in [Eqn (4.76)] of the book.
    % This matches the TNO solver, and is shown in Equations (6) and (7) of [3].
    %
    psi_dot_omega = wheel_omega*-sin_gamma;  % psi_dot due to wheel rotation and gamma [Eqn (4.E4) Page 177 - Book]
    eps_gamma = t.PECP1*(1 + t.PECP2*dFz); % [Eqn (4.90) Page 186 - Book] Camber reduction factor
    %assert(eps_gamma < 1); % eps_gamma == 1 will cause infinite stiffnesses
    % avoid division by zero with (1-eps_gamma)
    eps_gamma = min(eps_gamma, 1 - 1e-6);
    phi = phi_t + (1-eps_gamma) * psi_dot_omega * inv_Vcx_prime; % [Eqn (4.76) Page 184 - Book]

    % REVIEW: This limits turnslip forces for slow parking-lot maneuvers,
    % which we probably don't want...
    phi = phi * low_speed_scale;

    % normalized phi
    nphi = R0 * phi;

    % zeta3
    % cornering stiffness slope reduction factor
    zeta3 = cos_magic_formula_1(SQR(nphi), t.PKYP1); % [Eqn (4.79) Page 185 - Book]
  end
  %end
  %[phi, nphi, eps_gamma, zeta3] = compute_spin_slip(t, phi_t, sin_gamma, wheel_omega, inv_Vcx_prime, dFz, R0, want_turn_slip);
  
  %% Slip stiffnesses
  % These are required for relaxation lengths and force computations.
  %
  % This function computes 2 different cornering stiffnesses:
  %
  % Kya is the cornering stiffness: dFy0 / dalpha at alpha == 0 and kappa == 0 (i.e. pure Fy), but both gamma and phi possibly nonzero.
  % Kya_0 is dFy0 / dalpha at alpha, kappa, gamma and phi == 0. Used for combined Mz computation
  %
  %function [Kxk, Kya, Kya_0, inv_Kya] = calc_slip_stiffnesses(t, Fz, dFz, nFz_prime, gamma_star, dpi, zeta3, Fz0_prime)

  % Longitudinal slip stiffness Kxk

  % Kxk = BxCxDx = dFx0/dkappa_x at kappa_x = 0 (= CFk) (4.E15)
  Kxk = Fz * (t.PKX1 + t.PKX2*dFz) * exp(t.PKX3*dFz) * (1 + t.PPX1*dpi + t.PPX2*SQR(dpi)) * t.LKX;

  % Lateral slip stiffness: part of Fy computation
  Dya0 = t.PKY1 * Fz0_prime * (1 + t.PPY1*dpi);
  Dya  = Dya0 * (1 - t.PKY3*abs(gamma_star));

  Bya0 = 1 / ((t.PKY2 + 0                     ) * (1 + t.PPY2*dpi));
  Bya  = 1 / ((t.PKY2 + t.PKY5*SQR(gamma_star)) * (1 + t.PPY2*dpi));

  Cya0 = t.PKY4;
  Cya = Cya0;

  % Lateral slip stiffness Kya

  % Kya = dFy0 / dalpha at alpha == 0 but both gamma and phi possibly != 0
  % Kya_0 = dFy0 / dalpha at alpha == 0, gamma == 0, and phi == 0
  [Kya, Kya_0] = two_sin_magic_formula(nFz_prime,nFz_prime, Bya,Bya0, Cya,Cya0, Dya,Dya0, 0.0,0.0);
  %Kya = sin_magic_formula_3(nFz_prime, Bya, Cya, Dya);
  %Kya_0 = sin_magic_formula_3(nFz_prime, Bya0, Cya0, Dya0);
  Kya   = Kya * t.LKY;
  Kya_0 = Kya_0 * t.LKY;
  
  % Incorporate turn_slip factor zeta3
  % zeta3 is phi effect on cornering stiffness
  Kya = Kya * zeta3;

  inv_Kya = 1 / (Kya - 1e-6);
  %assert(Kya <= 0);
  %assert(Kya_0 <= 0);

  %end
  %[Kxk, Kya, Kya_0, inv_Kya] = calc_slip_stiffnesses(t, Fz, dFz, nFz_prime, gamma_star, dpi, zeta3, Fz0_prime);
  
  %% Relaxation lengths
  %function [sigma, CFx, CFy] = calc_sigma(t, Fz, dFz, nFz, gamma, Kxk, Kya, dpi, R0, Fz0_prime)
 
  PCX0 = t.LONGITUDINAL_STIFFNESS;
  CFx = PCX0 * (1 + t.PCFX1 * dFz + t.PCFX2 * SQR(dFz)) * (1 + t.PCFX3 * dpi); % (Paper Eq. 17)

  PCY0 = t.LATERAL_STIFFNESS;
  CFy = PCY0 * (1 + t.PCFY1 * dFz + t.PCFY2 * SQR(dFz)) * (1 + t.PCFY3 * dpi); % (Paper Eq. 18)

  if (t.model_version <= 520) %then

    % MF 5.2 relaxation length equations
    sigma_y_g = 1 - t.PKY3 * abs(gamma);
    Fz_over_PTY2 = Fz / (t.PTY2 * Fz0_prime);

    % relaxation lengths
    sigma = [
      (t.PTX1 + t.PTX2 * dFz) * exp(-t.PTX3*dFz) * t.LSGKP * R0*nFz % From the MF-Tyre equation manual
      sigma_y_g * sin_magic_formula_3(Fz_over_PTY2, 1.0, 2.0, t.PTY1) * R0*t.LFZO*t.LSGAL];

  else

    % relaxation lengths
    sigma = [Kxk/CFx; -Kya/CFy]; % (Paper Eq. 19, 20)
  end
  
  %end
  %[sigma, CFx, CFy] = calc_sigma(t, Fz, dFz, nFz, gamma, Kxk, Kya, dpi, R0, Fz0_prime, steady_state, transient_model);
  
%% Kappa and Alpha computations

  if (~isfinite(rolling_radius)) %then
    rolling_radius = t.compute_rolling_radius(unclamped_nFz, R_omega, CFz);
  end
  
  % Tire slip velocities at point S in W frame
  % (Vc(2) has been side-flipped)
  Vs = [-(Vc(1) - rolling_radius * wheel_omega); Vc(2)];
  
  % Nonlinear Transient Model 
  % [Section 7.2, 7.3 -- Book]
  % [Eqn (8.122-8.132) Page 396 - Book]
  
  %function [kappa, tan_alpha, u_dot, v_dot] = transient_model(t, u, v, Vs, Vc, inv_Vcx_prime, sigma)

  if ~isempty(states)
    % Semi-Non-Linear model
    % [Eqn (7.7, 7.9) Page 331 - Book]
    % Eq 90,91 of [3]
    % REVIEW: Camber and phi dynamics aren't modelled here (7.11, 7.12)
    % FIXME: contact_size(1)*phi term neglected (Eq 89 of [3])
    %
    % sigma(1) * u_dot + abs(Vc(1)) * u = sigma(1) * Vs(1)
    % sigma(2) * v_dot + abs(Vc(2)) * v = sigma(2) * Vs(2)
    %
    % kappa_prime = u / sigma(1) * sign_Vcx
    % tan_alpha_prime = v / sigma(2)
    % 
    assert(transient_model == 1);
    
    u = states(1);
    v = states(2) * side_flip;

    kappa = u / sigma(1);
    tan_alpha = v / sigma(2);
    
    u_dot = Vs(1) - abs(Vc(1)) * kappa;
    v_dot = Vs(2) - abs(Vc(1)) * tan_alpha;
  else
    
    kappa = Vs(1) * inv_Vcx_prime;
    tan_alpha = Vs(2) * inv_Vcx_prime;
    
    u_dot = 0;
    v_dot = 0;
  end
  
  %end
  %[kappa, tan_alpha, u_dot, v_dot] = transient_model(t, u, v, Vs, Vc, inv_Vcx_prime, sigma);
  alpha = atan(tan_alpha) * sign_Vcx;
  
  % Clamp values to ensure we don't extrapolate too far
  unclamped_kappa = kappa;
  kappa = CLAMP(unclamped_kappa, t.KPUMIN, t.KPUMAX);

  alpha = CLAMP(alpha, t.ALPMIN, t.ALPMAX);

  % Use of alpha_star definition (Large slip angles correction)
  alpha_star = alpha;
  if (use_alpha_gamma_star) %then
    alpha_star = tan_alpha*sign_Vcx; % [Eqn (4.E3) Page 177 - Book]
  end

  % Grip adjustment factors
  LMUX_star = road_grip * t.LMUX;
  LMUY_star = road_grip * t.LMUY;

  % Slippery surface with friction decaying with increasing (slip) speed
  % Grip slip speed effect
  if (t.LMUV ~= 0) %then
    norm_Vs = sqrt(SQR(Vs(1))+SQR(Vs(2)));
    muv_scale = 1/(1 + t.LMUV*norm_Vs*inv_V0);  % [Eqn (4.E7) Page 179 - Book]
    LMUX_star = LMUX_star * muv_scale;
    LMUY_star = LMUY_star * muv_scale;
  end

  % Digressive friction factor
  LMUX_prime = LMUX_star;
  LMUY_prime = LMUY_star;
  Amu = t.LAMU;
  if (Amu ~= 1) %then
    % On Page 179 of the book is suggested Amu=10, but after comparing the use
    % of the scaling factors against TNO, Amu = 1 was giving perfect match
    LMUX_prime = Amu * LMUX_star / (1 + (Amu-1) * LMUX_star); % [Eqn (4.E8) Page 179 - Book]
    LMUY_prime = Amu * LMUY_star / (1 + (Amu-1) * LMUY_star); % [Eqn (4.E8) Page 179 - Book]
  end
  
  %% Fx0 - Pure Longitudinal Force
  %function [Fx0, mu_x] = calc_Fx0(t, kappa, gamma, nphi, dFz, Fz, Kxk, LMUX_prime, LMUX_star, dpi, low_speed_scale, want_turn_slip)
  
  mu_x = max(0.0, (t.PDX1 + t.PDX2*dFz) * (1 + t.PPX3*dpi + t.PPX4*SQR(dpi)) * (1 - t.PDX3*SQR(gamma))*LMUX_star); % (4.E13)

  zeta1 = 1;
  if (want_turn_slip) %then
  
    % zeta1 - weight to diminish longitudinal force with increasing spin
    Bxp = t.PDXP1 * (1 + t.PDXP2*dFz) * cos_magic_formula_1(kappa, t.PDXP3); % [Eqn (4.106) Page 188 - Book]

    zeta1 = cos_magic_formula_1(nphi, Bxp); % [Eqn (4.105) Page 188 - Book]
    
  end
  
  % (alpha = 0) Page 179 of the book
  Cx = t.PCX1*t.LCX; % (> 0) (4.E11)

  Dx = mu_x * Fz * zeta1; % (> 0) (4.E12)
  
  %assert(Kxk >= 0);
  %assert(Cx*Dx >= 0);
  Bx = Kxk / (Cx*Dx + 1e-6); % (4.E16)

  SHx = (t.PHX1 + t.PHX2*dFz) * t.LHX; % (4.E17)
  SHx = SHx * low_speed_scale;

  kappa_x = kappa + SHx; % (4.E10)
  
  Ex = min(1.0, (t.PEX1 + t.PEX2*dFz + t.PEX3*SQR(dFz)) * (1 - t.PEX4*SGN(kappa_x)) * t.LEX); % (<=1) (4.E14)
  
  SVx = Fz * (t.PVX1 + t.PVX2*dFz) * t.LVX * LMUX_prime; % (4.E18)
  SVx = SVx * zeta1;
  SVx = SVx * low_speed_scale;
  
  Fx0 = sin_magic_formula(kappa_x, Bx, Cx, Dx, Ex) + SVx; % (4.E9)
  
  Fx0 = Fx0 * sign_Vcx;
  
  %end
  %[Fx0, mu_x] = calc_Fx0(t, kappa, gamma, nphi, dFz, Fz, Kxk, dpi, LMUX_prime, LMUX_star, low_speed_scale, want_turn_slip);

  %% Combined Fx
  %function Fx = calc_combined_Fx(t, Fx0, dFz, kappa, gamma_star, alpha_star, sign_Vcx)

  % (Combined Slip), Page 181 of the book
  Cxa = t.RCX1; % (4.E55)
  Exa = min(1.0, t.REX1 + t.REX2*dFz); % (<= 1) (4.E56)

  SHxa = t.RHX1; % (4.E57)
  Bxa = (t.RBX1 + t.RBX3*SQR(gamma_star)) * cos_magic_formula_1(kappa, t.RBX2) * t.LXAL; % (> 0) (4.E54)

  SHxap = SHxa + alpha_star*sign_Vcx;

  [Gxa0, Gxa] = two_cos_magic_formula(SHxa, SHxap, Bxa, Bxa, Cxa, Cxa, 1.0, 1.0, Exa, Exa);
  %Gxa0 = cos_magic_formula(SHxa,  Bxa, Cxa, 1.0, Exa);
  %Gxa  = cos_magic_formula(SHxap, Bxa, Cxa, 1.0, Exa);
  Gxa = Gxa / Gxa0;

  Fx = Gxa * Fx0; % (4.E50)

  %end
  %Fx = calc_combined_Fx(t, Fx0, dFz, kappa, gamma_star, alpha_star, sign_Vcx);
  
  %% Fy0 - Pure Lateral Force 
  % Fy0, mu_y -> used in combined Fy, output 
  % Fy0_0 -> used in combined Mz 
  % Gyk, zeta2 -> used by combined_Fy, combined Mz 
  % SHy, SVy -> used in pure Mz 
  % Kyp0, Kyg0 -> output stiffnesses 
  % ByCy -> used in pure Mz, combined Mzr, transient validation 
  % 
  % This function computes two different values of Fy0:
  %         Fy0 at (gamma, phi), used in Fy computations
  %         Fy0_0 at (gamma == 0, phi == 0), used in combined Mz computations
  %  
  %function [Fy0, Fy0_0, mu_y, Gyk, Gyk_0, zeta2, SHy, SVy, Kyp0, Kyg0, ByCy] = calc_Fy0( ... 
  %     t, Fz, Fz0_prime, kappa, sign_Vcx, tan_alpha, alpha_star, gamma_star, ... 
  %     Kya, Kya_0, inv_Kya, low_speed_scale, R0, ...
  %     nphi, eps_gamma, LMUY_star, LMUY_prime, dFz, dpi, want_turn_slip, use_alpha_gamma_star)

  % mu_y_0 is mu_y at gamma == 0.
  mu_y_0 = (t.PDY1 + t.PDY2 * dFz) * (1 + t.PPY3*dpi + t.PPY4*SQR(dpi)) * LMUY_star; % (4.E23) 
  mu_y   = mu_y_0 * (1 - t.PDY3*SQR(gamma_star)); % (4.E23) 
   
  mu_y_0 = max(0.0, mu_y_0);
  mu_y   = max(0.0, mu_y);

  % (zeta0 in book is equivalent to (1 - want_turn_slip) -- a flag to enable/disable turn slip calculations )

  SVyg = Fz*(t.PVY3 + t.PVY4*dFz)*gamma_star * t.LKYC * LMUY_prime; % (4.E28) at phi == 0

  zeta2 = 1;
  if (want_turn_slip) %then
    % zeta2
    Byp = t.PDYP1*(1 + t.PDYP2*dFz) * cos_magic_formula_1(tan_alpha, t.PDYP3); % [Eqn (4.78) Page 185 - Book]
    % peak side force reduction factor

    zeta2 = cos_magic_formula_1(abs(nphi) + t.PDYP4*sqrt(abs(nphi)), Byp); % [Eqn (4.77) Page 184 - Book]
    
    SVyg = SVyg * zeta2; % (4.E28)
  end

  % (kappa = 0) Page 179 of the book
  SHyk = t.RHY1 + t.RHY2*dFz; % (4.E65)

  %assert(Kya <= 0);

  Byk_common = cos_magic_formula_1(alpha_star - t.RBY3, t.RBY2) * t.LYKA; % (> 0) (4.E62)
  Byk = Byk_common * (t.RBY1 + t.RBY4 * SQR(gamma_star)); 

  Cyk = t.RCY1; % (4.E63)

  Eyk = min(1.0, t.REY1 + t.REY2*dFz); % (<=1) (4.E64)
  
  SHykk = SHyk + kappa;
  [Gyk0, Gyk] = two_cos_magic_formula(SHyk, SHykk, Byk, Byk, Cyk, Cyk, 1.0, 1.0, Eyk, Eyk);
  %Gyk0 = cos_magic_formula(SHyk,  Byk, Cyk, 1.0, Eyk);
  %Gyk  = cos_magic_formula(SHykk, Byk, Cyk, 1.0, Eyk);
  
  Gyk = Gyk / Gyk0;
  
  % IMPORTANT NOTE: The book says Gyk_0 should be computed with gamma == 0 and phi/phi_t == 0
  % in the computation of combined Mz (4.E74).  The computation above is computed with
  % phi == 0, but a gamma effect is included in Byk.  None of the other factors
  % (SHyk, Cyk, Eyk) involve either gamma or phi so are the same for both.
  %
  % If t.RBY4 happens to be zero (as in reference book parameters),
  % then Gyk_0 always == Gyk.

  Gyk_0 = Gyk;
  
  % The only difference between Gyk and Gyk_0 is the Byk_0 term.
  % A comparison with TNO gives better results when Gyk is used in both
  % Fy0 and Mz computations.
  %
% if (t.RBY4 ~= 0) %then
%   Byk_0 = Byk_common * (t.RBY1 + t.RBY4 * SQR(0.0));
%   [Gyk0_0,Gyk_0] = two_cos_magic_formula(SHyk, SHykk, Byk_0, Byk_0, Cyk, Cyk, 1.0, 1.0, Eyk, Eyk);
%   %Gyk0_0 = cos_magic_formula(SHyk,  Byk_0, Cyk, 1.0, Eyk); 
%   %Gyk_0  = cos_magic_formula(SHykk, Byk_0, Cyk, 1.0, Eyk); 
%   Gyk_0 = Gyk_0 / Gyk0_0;
% end

  SHy_0 = (t.PHY1 + t.PHY2*dFz) * t.LHY;
  
  if (t.model_version <= 520) %then

    % From MF-Tyre equation manual
    % (simplified, without pressure effect)
    Kya_52 = Fz * sin_magic_formula_3(Fz, 1 / (t.PKY2 * Fz0_prime), t.PKY4, t.PKY1) * t.LKYC;

    Kyg0 = (t.PHY3 * Kya_52 + Fz * (t.PVY3 + t.PVY4 * dFz)) * t.LKYC;
    SHy = SHy_0 + t.PHY3 * gamma_star * t.LKYC;
    %%asym
    
    Kyp0 = 0;

  else

    % Kyg0 = dFy0/dgamma at alpha = 0, gamma = 0 (= CFgamma) (4.E30)
    Kyg0 = Fz*(t.PKY6 + t.PKY7 * dFz) * (1 + t.PPY5*dpi) * t.LKYC;

    if (want_turn_slip) %then
      % zeta4 computations

      CHyp = t.PHYP1; % (>0) [Eqn (4.85) Page 186 - Book]
      DHyp = (t.PHYP2 + t.PHYP3*dFz) * sign_Vcx; % [Eqn (4.86) Page 186 - Book]
      EHyp = t.PHYP4; % (<=1) [Eqn (4.87) Page 186 - Book]

      % KyRp0 = spin force stiffness: dFy / d(1/R) at alpha,gamma,phi == 0.
      % Same as moment stiffness against gamma Kyg0, except with the (1-eps_gamma) factor applied to gamma.
      % (1-eps_gamma) guaranteed non-zero
      KyRp0 = Kyg0 / (1-eps_gamma); % Eqn (4.89) spin force stiffness 

      BHyp = KyRp0 / (CHyp*DHyp*Kya - 1e-6); %[Eqn (4.88) Page 186 - Book]
      %assert(Kya <= 0);
      %assert(CHyp*DHyp*Kya <= 0);

      SHyp = sin_magic_formula(nphi, BHyp, CHyp, DHyp, EHyp) * sign_Vcx;

      SHy  = SHy_0 + (SHyp - SVyg * inv_Kya); % [Eqn (4.84) Page 186 - Book -- zeta4 without '1 + ' term]

      Kyp0 = KyRp0 * R0;
    else

      SHy  = SHy_0 + (Kyg0 * gamma_star - SVyg) * inv_Kya; % (4.E27)
      %%asym
      Kyp0 = 0;

    end
  end
  SHy_0 = SHy_0 * low_speed_scale;
  SHy   = SHy   * low_speed_scale;

  SVy_0 = Fz * (t.PVY1 + t.PVY2*dFz) * t.LVY*LMUY_prime; % (4.E29)
  SVy   = SVy_0 * zeta2 + SVyg; % (4.E29)
  
  SVy_0 = SVy_0 * low_speed_scale;
  SVy   = SVy   * low_speed_scale;

  Cy    = t.PCY1*t.LCY; % (> 0) (4.E21)
  Cy_0  = Cy;

  Dy    = mu_y   * Fz * zeta2; % (4.E22)
  Dy_0  = mu_y_0 * Fz; % (4.E22)

  alpha_y   = alpha_star + SHy; % (4.E20)
  alpha_y_0 = alpha_star + SHy_0; % (4.E20)
  Ey_0  = min(1.0, (t.PEY1 + t.PEY2*dFz) * (1 + 0                      - (t.PEY3 + 0                ) * SGN(alpha_y_0)) * t.LEY); % (<=1)(4.E24)
  Ey    = min(1.0, (t.PEY1 + t.PEY2*dFz) * (1 + t.PEY5*SQR(gamma_star) - (t.PEY3 + t.PEY4*gamma_star) * SGN(alpha_y)   ) * t.LEY); % (<=1)(4.E24)
  %%asym

  %assert(Cy*Dy >= 0);
  %assert(Cy_0*Dy_0 >= 0);
  By    = Kya   / (Cy*Dy + 1e-6); % (4.E26)
  By_0  = Kya_0 / (Cy_0*Dy_0 + 1e-6); % (4.E26)

  if (t.model_version <= 520) %then
    Fy0   = sin_magic_formula(alpha_y, By, Cy, Dy, Ey) + SVy; % (4.E19)
    Fy0_0 = Fy0;
  else
    [Fy0,Fy0_0] = two_sin_magic_formula(alpha_y, alpha_y_0, By,By_0, Cy,Cy_0, Dy,Dy_0, Ey,Ey_0);
    %Fy0   = sin_magic_formula(alpha_y,   By,   Cy,   Dy,   Ey  ); % (4.E19)
    %Fy0_0 = sin_magic_formula(alpha_y_0, By_0, Cy_0, Dy_0, Ey_0)0; % (4.E19)
    Fy0   = Fy0 + SVy;
    Fy0_0 = Fy0_0 + SVy_0;
  end

  if (use_alpha_gamma_star) %then
    Fy0   = Fy0   * sign_Vcx;
    Fy0_0 = Fy0_0 * sign_Vcx;
  end
  
  ByCy = By*Cy;
  
  %end
  %[Fy0, Fy0_0, mu_y, Gyk, Gyk_0, zeta2, SHy, SVy, Kyp0, Kyg0, ByCy] = calc_Fy0( ... 
  %     t, Fz, Fz0_prime, kappa, sign_Vcx, tan_alpha, alpha_star, gamma_star, ... 
  %     Kya, Kya_0, inv_Kya, low_speed_scale, R0, ...
  %     nphi, eps_gamma, LMUY_star, LMUY_prime, dFz, dpi, want_turn_slip, use_alpha_gamma_star);
  
  %% Combined Fy
  %function Fy = calc_combined_Fy(t, Fy0, mu_y, Gyk, Fz, dFz, gamma_star, alpha_star, kappa, zeta2, low_speed_scale)

  BVyk = t.RVY6;
  CVyk = t.RVY5;
  
  DVyk = mu_y * Fz * (t.RVY1 + t.RVY2*dFz + t.RVY3*gamma_star);  % (4.E67)
  %%asym
  DVyk = DVyk * cos_magic_formula_1(alpha_star, t.RVY4);
  DVyk = DVyk * zeta2;

  SVyk = sin_magic_formula_3(kappa, BVyk, CVyk, DVyk) * t.LVYKA;  % (4.E66)
  SVyk = SVyk * low_speed_scale;

  Fy = Gyk * Fy0 + SVyk; % (4.E58)
  
  %end
  %Fy = calc_combined_Fy(t, Fy0, mu_y, Gyk, Fz, dFz, gamma_star, alpha_star, kappa, zeta2, low_speed_scale);
  
  %% Pure Mz aligning torque
  %function [alpha_r, alpha_t, Bt, Ct, Dt, Et] = calc_Mz0(...
  %    t, SHy, SVy, inv_Kya, Fz, sign_Vcx, dFz, alpha_star, gamma_star, nphi, ...
  %    dpi, LMUY_star, want_turn_slip, R0, inv_Fz0_prime)
  
  % (kappa = 0) Page 180 of the book
  SHt = t.QHZ1 + t.QHZ2*dFz + (t.QHZ3 + t.QHZ4*dFz)*gamma_star; % (4.E35)
  %%asym
  SHf = SHy + SVy * inv_Kya; % (4.E38)

  alpha_r = alpha_star + SHf; % = alphaf (4.E37)
  alpha_t = alpha_star + SHt; % (4.E34)

  % IMPORTANT NOTE: The following original equation (4.E43) was not matching the
  % TNO solver. The coefficient Dt affects the pneumatic trail (t) and the
  % self aligning torque (Mz).
  % It was observed that when negative inclination angles where used as
  % inputs, there was a discrepancy between the TNO solver and mfeval.
  % This difference comes from the term QDZ3, that in the original equation
  % is multiplied by abs(gamma_star), rather than gamma_star in the paper.
  % Equation (86) from the paper resulted in a perfect match with TNO.
  % Book definition:
  % Dt0 = Fz*(R0*inv_Fz0_prime) * (t.QDZ1 + t.QDZ2*dFz) * (1 - t.PPZ1*dpi) * t.LTR * signVcx; % (4.E42)
  % Dt  = Dt0*(1 + QDZ3*abs(gamma_star) + QDZ4*SQR(gamma_star)); % (4.E43)
  Dt = Fz*(R0*inv_Fz0_prime) * t.LTR * sign_Vcx ...
      * (t.QDZ1 + t.QDZ2*dFz) ...
      * (1 - t.PPZ1*dpi) ...
      * (1 + t.QDZ3*gamma_star + t.QDZ4*SQR(gamma_star)); % (Paper Eq. 86)
%%asym    
  zeta5 = 1;
  if (want_turn_slip) %then
    % zeta5: scale of Dt
    zeta5 = cos_magic_formula_1(abs(nphi), t.QDTP1); % [Eqn (4.91) Page 186 - Book]
    Dt = Dt * zeta5;
  end

  Ct = t.QCZ1; % (> 0) (4.E41)

  % IMPORTANT NOTE: The equation in the book for Bt uses the parameter QBZ6,
  % which is not found in standard .TIR files.  Also note that QBZ6 does not
  % appear in the full set of parameters on pages 190 and 615 of the book.
  % Instead equation (84) from the paper is used.
  % Book definition:
  %Bt= (t.QBZ1 + t.QBZ2*dFz + t.QBZ3*SQR(dFz)) * ...
  %    (1 +                     t.QBZ5*abs(gamma_star) + QBZ6*SQR(gamma_star) ) * t.LKY/t.LMUY_star; %(> 0)(4.E40)
  
  Bt = (t.QBZ1 + t.QBZ2*dFz + t.QBZ3*SQR(dFz)) * ...
       (1 + t.QBZ4*gamma_star + t.QBZ5*abs(gamma_star) ) * t.LKY/LMUY_star; %(> 0) (Paper Eq. 84)
%%asym     
  Et = (t.QEZ1 + t.QEZ2*dFz + t.QEZ3*SQR(dFz)) * ...
       (1 + (t.QEZ4 + t.QEZ5*gamma_star) * (2/pi) * atan(Bt*Ct*alpha_t)); % (<=1) (4.E44)
%%asym
  Et = min(1.0, Et);
  
  %end
  %[alpha_r, alpha_t, Bt, Ct, Dt, Et] = calc_Mz0(...
  %    t, SHy, SVy, inv_Kya, Fz, sign_Vcx, dFz, alpha_star, gamma_star, nphi, ...
  %    dpi, LMUY_star, want_turn_slip, R0, inv_Fz0_prime);
  
  
%% Mz - Combined Aligning Torque
  %function [alpha_r_eq, alpha_t_eq] = calc_alpha_eq(t, alpha_r, alpha_t, Kxk, inv_Kya, kappa)

  % (Combined slip), Page 182 of the book

  % IMPORTANT NOTE: The equations 4.E78 and 4.E77 are not used due to small
  % differences discovered at negative camber angles with the TNO solver.
  % Instead equations 80 and 81 from [3] are used.
  % alpha_r_eq = sqrt(SQR(alpha_r) + SQR(Kxk * inv_Kya * kappa) * SGN(alpha_r); % (4.E78)
  % alpha_t_eq = sqrt(SQR(alpha_t) + SQR(Kxk * inv_Kya * kappa) * SGN(alpha_t); % (4.E77)
  
  alpha_eq_shift = SQR(Kxk * inv_Kya * kappa);
  [tan_alpha_r,tan_alpha_t] = two_tan(alpha_r, alpha_t);
  theta_alpha_r = sqrt(SQR(tan_alpha_r) + alpha_eq_shift) * SGN(alpha_r); % (Paper Eq. 80)
  theta_alpha_t = sqrt(SQR(tan_alpha_t) + alpha_eq_shift) * SGN(alpha_t); % (Paper Eq. 81)
  [alpha_r_eq,alpha_t_eq] = two_atan(theta_alpha_r, theta_alpha_t);
  
  %end
  %[alpha_r_eq, alpha_t_eq] = calc_alpha_eq(t, alpha_r, alpha_t, Kxk, inv_Kya, kappa);
  
  %% Pneumatic Trail
  %function trail = calc_trail(t, alpha_t_eq, Bt, Ct, Dt, Et, cos_alpha_prime, Fy, gamma_star, dFz, R0, inv_Fz0, low_speed_scale)

  % Note that trail(1) is positive when the effective point of Fy application is *behind* the contact center.
  trail(1) = cos_magic_formula(alpha_t_eq, Bt, Ct, Dt, Et) * cos_alpha_prime;

  % IMPORTANT NOTE : the following equation is not written in any source, but "trail"
  % is multiplied by LFZO in the TNO solver. This has been empirically discovered.
  % REVIEW: Should trail(2) be scaled by LFZO too?
  trail(1) = trail(1) * t.LFZO;

  % IMPORTANT NOTE: trail(2) below ("s" in book and paper, Equation 4.E76) determines the
  % effect of Fx on Mz. The book uses Fz0_prime in the formulation, but the paper uses Fz0.
  % The equation (A56) from the paper has a better correlation with TNO.
  % trail(2) = R0*(t.SSZ1 + t.SSZ2*(Fy*inv_Fz0_prime) + (t.SSZ3 + t.SSZ4*dFz)*gamma_star)*t.LS; % (4.E76)
  trail(2) =   R0*(t.SSZ1 + t.SSZ2*(Fy*inv_Fz0)       + (t.SSZ3 + t.SSZ4*dFz)*gamma_star)*t.LS; % (Paper Eq. 82)
  %%asym

  % low_speed_scale is applied to Mz by dialing down the trail, rather than to Fx and Fy_prime.
  % (Mzr is scaled directly by low_speed_scale)
  trail = trail * low_speed_scale; % scale both x and y

  %end
  %trail = calc_trail(t, alpha_t_eq, Bt, Ct, Dt, Et, cos_alpha_prime, Fy, gamma_star, dFz, R0, inv_Fz0, low_speed_scale);
  
  %% Residual Mz
  %function Mzr = calc_Mzr(t, mu_y, LMUY_star, Fz, dFz, nFz_prime, ByCy, sign_Vcx, ...
  %      gamma, gamma_star, cos_alpha_prime, alpha_r_eq, ...
  %      Gyk, nphi, eps_gamma, zeta2, dpi, low_speed_scale, R0, want_turn_slip)

  Br = (t.QBZ9*t.LKY/LMUY_star + t.QBZ10 * ByCy); % preferred: qBz9 = 0 (4.E45)
  Dr = Fz*R0*LMUY_star * sign_Vcx;
  
  if (want_turn_slip) %then

    % zeta8 -> Drp
    Mzp_inf = t.QCRP1 * R0 * mu_y * Fz * sqrt(nFz_prime) * t.LMP; % [Eqn (4.95) Page 187 - Book]
    epsilon_mzp = 1e-6;
    Mzp_inf = max(Mzp_inf, epsilon_mzp); % Mzp_inf should be always > 0

    CDrp = t.QDRP1; % (>0) [Eqn (4.96) Page 187 - Book]
    DDrp = Mzp_inf / sin(pi/2 * CDrp); % [Eqn (4.94) Page 187 - Book]
    EDrp = t.QDRP2; % (<=1) [Eqn (4.97) Page 187 - Book]

    Kzgr0 = Fz*R0*t.LKZC*(t.QDZ8 + t.QDZ9*dFz + (t.QDZ10 + t.QDZ11*dFz) * abs(gamma)); %[Eqn (4.99) Page 187 - Book]

    % KzRp0 = Mz stiffness: dMz / (d(1/R), at alpha,gamma,phi == 0.
    % Same as Mz stiffness against gamma Kzgr0, except with 1/(1-eps_gamma) factor applied to gamma.
    
    % This is:
    % KzRp0 = Kzgr0 / (1-eps_gamma)
    % BDrp = KzRp0 / (CDrp*DDrp)
    
    BDrp = Kzgr0 / ( CDrp*DDrp * (1-eps_gamma) + 1e-6); % [Eqn (4.98) Page 187 - Book]
    %assert(CDrp*DDrp*(1-eps_gamma) >= 0);
    
    Drp = sin_magic_formula(nphi, BDrp, CDrp, DDrp, EDrp);
    
    % zeta6: scale of Br
    zeta6 = cos_magic_formula_1(abs(nphi), t.QBRP1); % [Eqn (4.102) Page 188 - Book]

    Br = Br * zeta6;

    Dr = Dr * (zeta2 * (t.QDZ6 + t.QDZ7*dFz) * t.LRES) + Drp; % (4.E47)
    
    % zeta7 == Cr == torque level at alpha = 90deg

    % IMPORTANT NOTE: Book uses R0*abs(phi_t) in Eq. (4.103) instead of R0*abs(phi) (nphi) as used in [3] and [4]
    
    % Mzp_90 = Mzp_inf * (2/pi) * atan(t.QCRP2 * abs(R0*phi_t)) * Gyk; % [Eqn (4.103) Page 188 - Book]
    Mzp_90   = Mzp_inf * (2/pi) * atan(t.QCRP2 * abs(nphi)) * Gyk; % (Page 25 of [4])

    epsilon_r = 1e-6;
    theta_Cr = CLAMP(Mzp_90 / (abs(Drp) + epsilon_r), -1.0, 1.0);
    Cr = (2/pi) * acos(theta_Cr); % [Eqn (4.104) Page 188 - Book]

  else

    %assert(zeta2 == 1);
    Dr = Dr * ( ...
      (t.QDZ6 + t.QDZ7*dFz) * t.LRES ...
       + gamma_star * ( (t.QDZ8 + t.QDZ9*dFz) * (1 + t.PPZ2*dpi) + (t.QDZ10 + t.QDZ11*dFz) * abs(gamma_star) ) * t.LKZC ...
      ); % (4.E47)
    %%asym
    
    Cr = 1;
  end
  
  Mzr = cos_magic_formula_3(alpha_r_eq, Br, Cr, Dr) * cos_alpha_prime; % (4.E36, 4.E75)
  
  Mzr = Mzr * low_speed_scale;
  
  %end
  %Mzr = calc_Mzr(t, mu_y, LMUY_star, Fz, dFz, nFz_prime, ByCy, sign_Vcx, ...
  %      gamma, gamma_star, cos_alpha_prime, alpha_r_eq, ...
  %      Gyk, nphi, eps_gamma, zeta2, dpi, low_speed_scale, R0, want_turn_slip);
  
  %% Combined Mz
  %function Mz = calc_combined_Mz(t, Fy0_0, Gyk_0, Fx, trail, Mzr)
  
  Fy_prime = Gyk_0 * Fy0_0; % (4.E74)
  
  Mzx =  trail(2) * Fx;  % (4.E71)
  Mzy = -trail(1) * Fy_prime;
  
  Mz = Mzx + Mzy + Mzr;  % (4.E31)
  %end
  %Mz = calc_combined_Mz(t, Fy0_0, Gyk_0, Fx, trail, Mzr);
  
  %% Mx - Overturning Couple
  % (Also see Section 4.3.5), Page 182 of the book and Page 12 of the paper

  % IMPORTANT NOTE: Equation (4.E69) in the book is not used; instead
  % equation (49) of [2] is used, which matches the results of the official TNO solver.
  %
  % Mx = R0*Fz*(QSX1*LVMX - QSX2*gamma*(1 + PPMX1*dpi) + QSX3*(nFy)...  % (4.E69)
  %     + QSX4*cos(QSX5*atan((QSX6*(nFz))^2))*sin(QSX7*gamma + QSX8*atan...
  %     (QSX9*(nFy))) + QSX10*atan(QSX11*(nFz))*gamma)*LMX; %(4.E69)
  %
  % IMPORTANT NOTE: In the book, the atan() term involving QSX6 is missing
  % a set of parentheses:  It is written (atan(QSX6*nFz))^2 rather than
  % atan( (QSX6)*nFz)^2 ) as shown in Eq (49) of [2].
  %
  % Book:   cos(QSX5 * atan(QSX6 * Fz / Fz0)) ^ 2)
  % Paper:  cos(QSX5 * atan((QSX6 * Fz / Fz0) ^ 2))   % odd that QSX6 is inside parentheses
  % I think the correct definition is:
  % Likely: cos(QSX5 * atan(QSX6 * (Fz / Fz0) ^ 2))   % This seems most likely.  Misprint in book is ^2 should be inside, not outside paren.
  %     This is a typical cos_magic_formula form with QSX6 = B, as documented in the definition of QSX6 above
  %     Squaring the atan as printed in the book is also not a typical magic_formula formulation.
  %
  % NOTE: The complexity of the Mx formulas is presumably to handle the complicated
  %       case of high-camber motorcycle tires.  For typical low-camber tire conditions,
  %       the complexity is probably not worth the considerable computation cost.
  % 
  % QSX1,QSX2,QSX3 are the MF 5.2 model parameters.
  % QSX4-QSX11 were added for MF 6.0. QSX4-QSX8 are quite expensive.
  % QSX12-14 were added for MF 6.1.  QSX11 costs an atan().
  %
  % Mx is a rich source of left-right asymmetries involving gamma and Fy.
  
  %function Mx = calc_Mx(t, Fy, Fz, gamma, dpi, nFz, inv_Fz0, R0)
  
  % [2] Eq. 49
  nFy = Fy*inv_Fz0;
    
  % MF 5.2 equations with pressure correction:
  Mx_Fz = t.QSX1*t.LVMX - t.QSX2 * gamma * (1 + t.PPMX1 * dpi) + t.QSX3 * nFy;

  % QSX4-QSX9
  if (t.QSX4 ~= 0) %then
    [a6,a9]  = two_atan( SQR(t.QSX6 * nFz), t.QSX9 * nFy);
%     %a6 = SQR(atan(t.QSX6 * nFz));  % Book 4.E69
%     %a6 = atan(t.QSX6 * SQR(nFz));  % Magic formula, QSX6 = B
%     a6 = atan(SQR(t.QSX6 * nFz));   % Eq. 49 of [2]
%     a9 = atan(t.QSX9 * nFy);
    c5  = cos( t .QSX5 * a6);
    s78 = sin( t.QSX7 * gamma + t.QSX8 * a9 );
    
    Mx_Fz = Mx_Fz + t.QSX4 * c5 * s78;
  end
  if (t.QSX10 ~= 0) %then
    Mx_Fz = Mx_Fz + t.QSX10 * atan(t.QSX11 * nFz) * gamma;
  end
  Mx_Fz = Mx_Fz - t.QSX12 * gamma * abs(gamma);

  Mx_Fy = t.QSX13 + t.QSX14*abs(gamma);

  Mx = R0 * t.LMX * (Fz * Mx_Fz + Fy * Mx_Fy);
  
  Mx = Mx * low_load_scale;
  %end
  %Mx = calc_Mx(t, Fy, Fz, gamma, dpi, nFz, inv_Fz0, R0);
  
  %% My - Rolling Resistance Moment
  % (See Eqns (9.236, 9.230, 9.231)), Page 182 of the book
  %
  % IMPORTANT NOTE: Equation (4.E70) is not used because a) the first term should be Fz0 instead of Fz,
  % and b) the equation as printed gives positive My values instead of the expected negative My.
  % Instead, Eq. 50 of [2] is used.
  %
  % My = Fz*R0*(QSY1 + QSY2*(Fx*inv_Fz0) + QSY3*abs(Vcx/V0) + QSY4*(Vcx/V0)^4 ...  (4.E70)
  %     + (QSY5 + QSY6*(Fz*inv_Fz0))*SQR(gamma)) * ((Fz*inv_Fz0)^QSY7 * (p/pi0)^QSY8)*LMY;

  %function My = calc_My(t, Fx, gamma, nFz, sign_Vcx, Vc, dpi, inv_Fz0, inv_V0, R0, Fz0)
    nFx = Fx * inv_Fz0;
    nVcx = Vc(1) * inv_V0;
    % MF 5.2: These equations work with MF5.2 when QSY5,6,7,8 are all 0
    My = -R0 * Fz0 * t.LMY * ( ... % Paper [2] Eq. 50
      t.QSY1 ...
      + t.QSY2 * nFx ...
      + t.QSY3 * abs(nVcx) ...
      + t.QSY4 * SQR(SQR(nVcx)) ...
      + (t.QSY5 + t.QSY6 * unclamped_nFz) * SQR(unclamped_gamma) ...
      );
    if (t.QSY7 ~= 0) %then
      My = My * POW(unclamped_nFz, t.QSY7);
    end
    if (t.QSY8 ~= 0) %then
      My = My * POW(dpi+1, t.QSY8);
    end
    My = My * sign_Vcx;
    
    My = My * low_load_scale;
    
  %end
  %My = calc_My(t, Fx, gamma, nFz, sign_Vcx, Vc, dpi, inv_Fz0, inv_V0, R0, Fz0);
  
    % Low-speed limit checking of My disabled for better comparison with MFeval
%     if (0) %then
%       if (low_speed_scale < 1) %then
%         s = abs(Vc(1))/t.VXLOW;
%         % sharp tapering to 0: d/ds==1 at s==0, d/ds=0 at s==1.
%         My = My * SHARPSTEP(s);
%       end
%     else  
%       % From Solver.m
%       % Low-speed reduction of rolling resistance
%       % (Empirically discovered)
%       kappa_hi = t.VXLOW/abs(Vcx_prime) - 1;
%       kappa_lo = -1 - t.VXLOW - kappa_hi;
%       % Low speed model (Empirically discovered)
%       if (kappa_lo <= kappa && kappa <= kappa_hi) %then
%           s = LERP(kappa, -1.0, kappa_hi, 0.0, 1.0);
%           My = My * sin(s*pi/2);
%       end
%     end            

  %% Compute loaded radius from computed loads
  if (steady_state) %then
    
    % Compute loaded radius, stiffness, and rho from vertical force
    [loaded_radius, stiffness, rho] = t.loaded_radius_from_load( ...
        unclamped_Fz, R_omega, n_omega, Fx, Fy, gamma, CFz, dpi);

    % Compute contact patch location
    radial = -loaded_radius * R(:,3);
    tire_contact = [radial(1:2);0];
    %center = [0;0;-radial(3)]; % contact - radial
  end

  %% Validate transient derivatives
  if (~isempty(states)) %then
    assert(transient_model == 1);
    
    alpha_eq = alpha_r_eq;
    alpha_sl = 3.0/ByCy; % == 3*Dy * inv_Kya;

    % at low speed, limit transient kappa and tan_alpha
    % by setting derivatives to zero if needed (Eq. 7.25)
    if (abs(Vc(1)) < t.VXLOW && abs(alpha_eq) > alpha_sl) %then
      if (states * kappa < 0) %then
        u_dot = 0;
      end
      if (v_dot * tan_alpha < 0) %then
        v_dot = 0;
      end
    end
    
    states_dot = [u_dot; v_dot * side_flip];
  else
    states_dot = [];
  end
  
  %% Compute contact patch size using vertical tire stiffness
  
  %function contact_size = calc_contact_patch(t, unclamped_Fz, unclamped_nFz, dpi, R0)
  if (t.model_version <= 520) %then

    % Note that we use t.Q_A1_ and t.Q_A2_ here,
    % computed in on_parameters_changed().

    % Times two to reflect full length and width
    contact_size = [
      2 * R0 * (t.Q_A2_*unclamped_nFz + t.Q_A1_ * sqrt(unclamped_nFz)) % From the MF-Tyre equation manual
      t.WIDTH]; % From the MF-Tyre equation manual

  else
    
    % REVIEW: Appendix 3 of [1], Eq (A3.7, A3.8) says that
    % the contact patch dimensions are best modeled as a function
    % of rho rather than the estimate unclamped_Fz/CFz which
    % uses the constant stiffness CFz.
    %
    % Which to use depends on the formulas used to perform the fitting,
    % though the results obtained either way should be quite close.
    %
    ss = unclamped_Fz / CFz;
    %ss = rho;
    ss = max(0.0, ss);  % in case of no_wheel_lift
    % Times two to reflect full length and width
    contact_size = [
      2 * R0      * (t.Q_RA2*ss + t.Q_RA1 * sqrt(ss))           % Paper Eq. 9
      2 * t.WIDTH * (t.Q_RB2*ss + t.Q_RB1 * POW(ss, 1.0/3.0))]; % Paper Eq. 10
  end
  
  %end
  %contact_size = calc_contact_patch(t, unclamped_Fz, unclamped_nFz, dpi, R0);

%% Store output

  if (low_speed_scale < 1) %then
    Fx = Fx * low_speed_scale;
    Fy = Fy * low_speed_scale;
    Mx = Mx * low_speed_scale;
    My = My * low_speed_scale;
    %% NO: Mz = Mz * low_speed_scale;
% NOTE: low_speed_scale has already been applied to Mz and trail
% (we want to be sure we don't apply low_speed_scale twice)
  end

  if (side_flip == -1) %then
    Fy = -Fy;
    Mz = -Mz;
    Mx = -Mx;

    alpha = -alpha;
    gamma = -gamma;
    
    phi_t = -phi_t;

    trail(2) = -trail(2);
    
    Vs(2) = -Vs(2);
    tire_contact_vel(2) = -tire_contact_vel(2);
  end

  output.right_tire = right_tire;

  output.R = R;
  output.W = W;

  output.kappa = kappa;
  output.alpha = alpha;
  output.gamma = gamma;
  output.path_curvature = phi_t;
  output.phi = phi;

  output.rolling_radius = rolling_radius;
  output.loaded_radius = loaded_radius;
  output.rho = rho;
  %output.pressure = pressure;

  output.tire_contact = tire_contact;
  output.tire_contact_vel = tire_contact_vel;
  output.tire_slip_vel = Vs;

  output.wheel_omega = wheel_omega;

  output.contact_size = contact_size;
  output.vertical_stiffness = stiffness;

  output.force = [Fx; Fy; unclamped_Fz];
  output.moment = [Mx; My; Mz];

  if ~isempty(algebraic_loops)
    algebraic_loops = [Fx; Fy];
    if use_adapted_SAE
      algebraic_loops(2) = -algebraic_loops(2);
    end
  end

  % dFx/dkappa
  % dFy/dalpha
  % dFy/dgamma at gamma == 0
  % dMz/dphi at gamma == 0
  output.slip_stiffness = [Kxk; Kya; Kyg0; Kyp0];
  output.relaxation = sigma;
  output.trail = trail;
  output.grip = [mu_x; mu_y];

  % Map output values if using "adapted SAE" units
  if (use_adapted_SAE) %then
    output.force(2) = -output.force(2);
    output.moment(3) = -output.moment(3);
    output.phi = -output.phi;
    output.path_curvature = -output.path_curvature;
    output.tire_slip_vel(2) = -output.tire_slip_vel(2);
    output.tire_contact_vel(2) = -output.tire_contact_vel(2);
    output.W(:,2) = -output.W(:,2);
    output.R(:,2) = -output.R(:,2);
  end
end

%% Fixed loaded radius eval()
% Steady

% solve for Fz at given loaded_radius
function output = fixed_loaded_radius_eval(t, input, output, options)

  fixed_radius = input.loaded_radius;
  input.loaded_radius = nan;
  
  options.no_wheel_lift = true;
  
  function try_delta_radius = try_fun(try_Fz)
    input.Fz = try_Fz;
    output = t.compute_steady(input, options);
    try_delta_radius = output.loaded_radius - fixed_radius;
  end
  
  ytol = t.UNLOADED_RADIUS * 1e-6;
  
  % use input Fz as starting point
  initial_Fz = max(50, input.Fz); % ensure non-zero
  range = 0.10;
  
  % if not provided, use nominal Fz and use wider starting range
  if (~isfinite(initial_Fz)) %then
    initial_Fz = t.FNOMIN;
    range = 0.50;
  end
  
  Fz0 = initial_Fz * (1 - range);
  Fz1 = initial_Fz * (1 + range);
  
  solve_secant(@try_fun, Fz0, Fz1, ytol);
  
  % last call to try_fun will set up output for return
end

%% Compute loaded radius from load
function [loaded_radius, stiffness, rho] = loaded_radius_from_load(t, Fz, R_omega, n_omega, Fx, Fy, gamma, CFz, dpi)
  % REVIEW: Simple quadratic solution won't work when bottoming is enabled.
  % In that case, might as well just do iterative solve.
  if (t.model_version < 611) %then

    % For model parameters < MF 6.2, we can calculate loaded radius directly from the rho quadratic
    [KR1, KR2, ~] = t.compute_rho_quadratic(R_omega, n_omega, Fx, Fy, gamma, dpi, false, 0);
    
    % Solve quadratic KR2*SQR(rho) + KR1*rho = Fz for rho:
    if (KR2 == 0) %then
      rho = Fz / KR1;
    else
      %assert(4*KR2*Fz + SQR(KR1) >= 0);
      rho = (sqrt(SQR(KR1) + 4*KR2*Fz) - KR1) / (2*KR2);
    end
    
    loaded_radius = R_omega - rho;

  else
    
    % For MF 6.2 and above we must find an iterative solution
    % because of the coupling of loaded_radius/rho, Fy, gamma
    % and the quadratic coefficients.
    
    lr0 = R_omega - Fz/CFz; % an initial guess based on nominal stiffness
    Fz_tol = t.FNOMIN * 1e-6;
    loaded_radius = solve_secant(@delta_fz_func, lr0 * 0.97, lr0 * 1.03, Fz_tol);
    % rho and q are set up by final call to delta_fz()
  end
  
  % Compute stiffness before return
  stiffness = KR1 + 2*KR2*rho;
  return;
  
  function delta_fz = delta_fz_func(try_radius)
    no_wheel_lift = true; % in case of negative Fz, e.g. fixed loaded radius computations
    [KR1, KR2, rho] = t.compute_rho_quadratic(R_omega, n_omega, Fx, Fy, gamma, dpi, no_wheel_lift, try_radius);
    
    try_Fz = (KR1*rho + KR2*SQR(rho));
    delta_fz = try_Fz - Fz;
  end
  
end

%% Compute coefficients of a quadratic in rho:
%   Fz = KR1*rho + KR2*rho^2
%
function [KR1, KR2, rho] = compute_rho_quadratic(t, R_omega, n_omega, Fx, Fy, gamma, dpi, no_wheel_lift, loaded_radius)

  R0 = t.UNLOADED_RADIUS;
  inv_R0 = 1/R0;
  Fz0 = t.FNOMIN;
  inv_Fz0 = 1/Fz0;
  
  % NOTE: rho is used below only if 6.2 model is in use
  rho = R_omega - loaded_radius;
  
  % If new 6.2 loaded radius model is in use...
  rho_zg = 0; % Radius adjustment due to camber
  if (t.Q_CAM3 ~= 0) %then
    
    % Improved loaded radius model from paper [4]: http://www.mate.tue.nl/mate/pdfs/11495.pdf
    % Section 5 describes the loaded radius model improvements.
    % Model supported in MF-Tire 6.2 and upwards.
    
    % loaded_radius required for 6.2 model
    if (gamma ~= 0) %then
      
      % reference tread width: similar to E.T.R.T.O definition
      rtw = t.WIDTH * (1.075 - 0.5 * t.ASPECT_RATIO);
      % half rtw tan gamma
      abs_gamma = abs(gamma);
      h = (rtw/2) * tan(abs_gamma);

      % limit magnitude of Q_CAM3 term
      rho_h = max(rho, -h);
      
      % Solver.m variant: hard transition at rho == 0
      num = SQR( (t.Q_CAM1*loaded_radius + t.Q_CAM2*SQR(loaded_radius)) * abs_gamma);
      den = SQR( (t.Q_CAM1*R_omega       + t.Q_CAM2*SQR(R_omega)      ) * abs_gamma);
      off = (t.Q_CAM3*abs_gamma * rho_h);
      rho_zg = (num * h) / (den * 4) - off;
      
%       % Model from [4], with soft transition at rho == 0
%       if (rho > 0) %then
%         % loaded_radius <= R_omega
%         num = SQR( (t.Q_CAM1*loaded_radius + t.Q_CAM2*SQR(loaded_radius)) * abs_gamma) - (t.Q_CAM3*abs_gamma * loaded_radius);
%         den = SQR( (t.Q_CAM1*R_omega       + t.Q_CAM2*SQR(R_omega)      ) * abs_gamma) - (t.Q_CAM3*abs_gamma * R_omega);
%         %assert(den > 0);
%         rho_zg = (num * h) / (den * 4);
%       else % rho <= 0
%         % loaded_radius > R_omega && loaded_radius < R_omega + h
%         rho_zg = SQR(h + rho_h) / (h * 4);
%       end
    end
  end
  
  fy_effect = 0;
  if (t.Q_FCY ~= 0) %then
    
    % Asymmetric effect on stiffness for combinations of camber and lateral force
    nlr = loaded_radius/R0;
    SFyg = (t.Q_FYS1 + t.Q_FYS2*nlr + t.Q_FYS3*SQR(nlr)) * gamma;
    %%asym, but follows sign of Fy
    
    Q_FCY_ = t.Q_FCY;
    if (t.Q_FCY2 ~= 0) %then
      % IMPORTANT NOTE: In the paper, Q_FCY2 is the exponent of
      % rho/R0, but in some MF 6.2 data (origin unknown)
      % the value of Q_FCY is NEGATIVE (-.4751), which results in
      % huge values at small rho (and inf at rho == 0).
      % Applying this exponent to loaded_radius/R0
      % is consistent with SFyg formulation, and gives reasonable results.
      % 
      %Q_FCY_ = Q_FCY_ * POW(rho*inv_R0, t.Q_FCY2); % Eq. 5.26 of [4]
      Q_FCY_ = Q_FCY_ * POW(nlr, t.Q_FCY2); % Modified Eq 5.26 of [4] 
    end
    fy_effect = SQR(   Q_FCY_ * inv_Fz0 * (Fy - SFyg) );
  end
  
  rho = rho + rho_zg;
  if (~no_wheel_lift) %then
    rho = max(0.0, rho);  % Eq. 5.24 of [4]
  end

  speed_effect    = t.Q_V2 * abs(n_omega);
  fx_effect       = SQR( t.Q_FCX  * inv_Fz0 * Fx);
  pressure_effect = (1 + t.PFZ1*dpi);
  
  % IMPORTANT NOTE: In Eq. 5.26 of [4], there is a missing set of parentheses
  % around the first term of Fcorr.  The formulation in Eq. 5.12 earlier
  % in the paper is correct.
  Fcorr = (1 + speed_effect - fx_effect - fy_effect) * pressure_effect; % [Eq 5.12 of [4]]
  %assert(Fcorr >= .5);
  
  Fcorr = Fcorr * Fz0 * inv_R0;
  
  % Coefficients of quadratic in rho
  % derived from Fz = Fcorr*Fz0 * (Q_FZ1*(rho/R0) + t.Q_FZ2*SQR(rho/R0));
  % Note that t.Q_FZ1_ is precomputed in on_params_changed();
  KR1 = Fcorr * (t.Q_FZ1_ + (t.Q_CAM + t.Q_FZ3) * SQR(gamma));  % include obsolete Q_FZ3 for back compatibility
  KR2 = Fcorr * inv_R0 * t.Q_FZ2;
end

%% Function to determine whether there is an algebraic loop
% involving Fx or Fy, and the force() indices affected.
%
function coupled = has_algebraic_loop(t)
  if (t.Q_FCX ~= 0 && t.Q_FCY ~= 0)
    coupled = [1;2];
  elseif (t.Q_FCX ~= 0) %then
    coupled = 1;
  elseif (t.Q_FCY ~= 0) %then
    coupled = 2;
  else
    coupled = [];
  end
end

%% Function to close algebraic loops
% Assumes both fx and fy are coupled.
% Usually one Newton-Raphson step is enough: a total of
% 3 extra evaluations.
function [input, output] = close_algebraic_loop(t, input, output, options, f_tol)
  if nargin < 5, f_tol = 30; end
  
  sel = t.has_algebraic_loop();
  if (isempty(sel)) %then
    return;
  end
  
  y = input.force(sel) - output.force(sel);
  if (y'*y < f_tol*f_tol) %then
    return;
  end
  
  x = output.force(sel);
  while true
    y = try_close(x);
    if (y'*y < f_tol*f_tol) %then
      return;
    end
    
    % Newton-Raphson
    J = one_sided_jacobian(@try_close, x, y);
    
    dx = J\y;
    %assert(norm(dx) > 0);
    
    x = x - dx;
  end
  
  function y = try_close(x)
    input.force(sel) = x;
    output = t.eval(input, options);
    y = input.force(sel) - output.force(sel);
  end
end

end % methods(public)

methods(Access = protected)

% Called after properties have changed (e.g. after load)
% (should be private, but used internally)
function on_params_changed(tire)
  
  % Don't do anything if parameters aren't valid
  if (~tire.is_valid) %then
    return;
  end
  
  tire.right_tire = (tire.TYRESIDE ~= 0);
  
  % Set model_version based on FITTYP
  fittyp = tire.FITTYP;
  % Some special cases 
  model = 0;
  switch (fittyp) %then
    case { 5, 51 } % 5 = MF-Tire 5.1, 51 = MF-MCTire 1.0
      model = 510;
    case { 6, 21, 52 } % 52 = MF-MCTire 1.1
      model = 520;
    otherwise
      if (fittyp >= 60 && fittyp <= 69)
        model = fittyp * 10;
      end
  end
  if (model < 520) %then
    fprintf('FITTYP = %d Model NOT SUPPORTED, assuming 5.2\n', fittyp);
    model = 520;
  end
  
  tire.model_version = model;

  % Compute the constant Q_FZ1 loaded radius model term.
  % This is computed from CFz0 = t.VERTICAL_STIFFNESS, which is defined as the
  % vertical tire stiffness at Fz = FNOMIN, pressure = NOMPRES,
  % and the parameters camber, Vcx, Fx, Fy all zero.
  %
  % We set the derivative of the loaded radius formula equal to CFz0,
  % plug in the conditions above, and solve for the unknown Q_FZ1.

  R0 = tire.UNLOADED_RADIUS;
  Fz0 = tire.FNOMIN;
  CFz0 = tire.VERTICAL_STIFFNESS;
  tire.Q_FZ1_ = sqrt(SQR(CFz0*R0/Fz0) - 4*tire.Q_FZ2);
  
  % If needed, compute default values for MF 5.2 contact length model params.
  if (model <= 520) %then
      tire.Q_A1_ = tire.Q_A1; % MF 5.2 Square root load term in contact length
      tire.Q_A2_ = tire.Q_A2; % MF 5.2 Linear load term in contact length
      if (tire.Q_A1_ == 0 && tire.Q_A2_ == 0) %then
          % Set default values (Empirically discovered)
          y = log10(R0*CFz0/Fz0);
          tire.Q_A1_ = y*(y*(y*-0.0388 + 0.2509) + -0.6283) + 0.6279;
          tire.Q_A2_ = 1.693 * SQR(tire.Q_A1_);
      end
  end
end

end % methods(Access = protected)

methods(Access = private)

%% Calculate effective rolling radius
function rolling_radius = compute_rolling_radius(t, unclamped_nFz, R_omega, CFz)
  Fz0 = t.FNOMIN;
  ndz = (Fz0/CFz);

  rr_nFz = max(0.0, unclamped_nFz); % No Fz effect if wheel is off the ground
  rolling_radius = R_omega - ndz * (t.DREFF * atan(t.BREFF*rr_nFz) + t.FREFF*rr_nFz); % [Eqn (7) Paper]
end

%% Compute effective rolling radius and an estimated omega
function [rolling_radius, wheel_omega] = compute_rolling_radius_and_wheel_omega(t, unclamped_Fz, kappa, Vcx, CFz)
  
  % Solve for a rolling radius x equal to the rolling radius
  % computed from Vcx, kappa and x.
  
  R0 = t.UNLOADED_RADIUS;
  V0 = t.LONGVL;
  inv_V0 = 1/V0;
  
  % Factor out values that don't change during iteration
  % (see compute_rolling_radius() and R_omega definition above)
 
  KRR_0 = (kappa + 1) * Vcx;
  KRR_1 = R0*t.Q_RE0;
  KRR_2 = R0*t.Q_V1*SQR(R0*inv_V0);
  
  Fz0 = t.FNOMIN;
  ndz = (Fz0/CFz);
  unclamped_nFz = unclamped_Fz / Fz0;
  
  rr_nFz = max(0.0, unclamped_nFz); % no Fz effect if wheel is off the ground
  KRR_3 = ndz * (t.DREFF * atan(t.BREFF*rr_nFz) + t.FREFF*rr_nFz); % [Eqn (7) Paper]
  
  % NOTE: This solve_secant()-based implementation is noticeably slower
  % than a simple inline iteration loop in standard Matlab, even
  % though it converges much faster.  This is due to the extra
  % nested-function call overhead.
  % In C++ (and presumably compiled Matlab) solve_secant()
  % has very little overhead.
  
  R_tol = KRR_1 * 1e-6;
  solve_secant(@delta_rolling_radius, KRR_1*0.95, KRR_1, R_tol);

  % Last call to delta_rolling_radius() will set
  % rolling_radius and wheel_omega for function return
  
  function delta_radius = delta_rolling_radius(try_radius)
    wheel_omega = KRR_0 / try_radius; % [Eqn (2.5) Page 65 - Book]
    R_omega = KRR_1 + KRR_2 * SQR(wheel_omega);    
    rolling_radius = R_omega - KRR_3;
    delta_radius = try_radius - rolling_radius;
  end
  
end

% Compute SWIFT model carcass stiffnesses and damping
% (not currently used)
function swift = compute_swift_params(t)
  % natural frequencies in rads/sec
  wn_xz  = t.FREQ_LONG*2*pi;
  wn_y   = t.FREQ_LAT*2*pi;
  wn_wup = t.FREQ_WINDUP*2*pi;
  wn_yaw = t.FREQ_YAW*2*pi;
  
  % carcass stiffnesses
  swift.kb_xz  = SQR(wn_xz)  * t.BELT_MASS; % N/m
  swift.kb_y   = SQR(wn_y)   * t.BELT_MASS; % N/m
  swift.kb_wup = SQR(wn_wup) * t.BELT_IYY;  % Nm/rad
  swift.kb_yaw = SQR(wn_yaw) * t.BELT_IZZ;  % Nm/rad
  
  % carcass damping
  swift.cb_xz  = t.DAMP_LONG   / (wn_xz   * 2*t.BELT_MASS); % N/(m/s)
  swift.cb_y   = t.DAMP_LAT    / (wn_y    * 2*t.BELT_MASS); % N/(m/s)
  swift.cb_wup = t.DAMP_WINDUP / (wn_wup * 2*t.BELT_IYY);  % Nm/(rad/s)
  swift.cb_yaw = t.DAMP_YAW    / (wn_yaw  * 2*t.BELT_IZZ);  % Nm/(rad/s)
end

end % end methods

properties
  
%[MODEL]
FITTYP  = 61            ;% 1 Magic Formula version number
PLUS    = 0             ;% 2 Moment representation 0=ground frame 1=wheel frame (not supported)
TYRESIDE = 0            ;% 3 Position of tire during measurements
LONGVL  = 16.7          ;% 4 Reference speed
VXLOW   = 1             ;% 5 Lower boundary velocity in slip calculation
ROAD_INCREMENT = 0.01   ;% 6 Increment in road sampling (Swift)
ROAD_DIRECTION = 1      ;% 7 Direction of travelled distance (Swift)
PROPERTY_FILE_FORMAT = 0 ;% 8 Tire model selection (Adams only)
USE_MODE = 0            ;% 9 Tire use mode switch (Adams only)
HMAX_LOCAL = 0          ;% 10 Local integration time step (Adams only)
TIME_SWITCH_INTEG = 0   ;% 11 Time when local integrator is activated (Adams only)
USER_SUB_ID = 0         ;% 12 Unknown (Adams only)
N_TIRE_STATES = 0       ;% 13 Unknown (Adams only)

%[DIMENSION]
UNLOADED_RADIUS = 0.3135 ;% 14 Free tire radius
WIDTH   = 0.205         ;% 15 Nominal section width of the tire
RIM_RADIUS = 0.1905     ;% 16 Nominal rim radius
RIM_WIDTH = 0.1651      ;% 17 Rim width
ASPECT_RATIO = 0.6      ;% 18 Nominal aspect ratio

%[OPERATING_CONDITIONS]
INFLPRES = 220000       ;% 19 Tire inflation pressure
NOMPRES = 220000        ;% 20 Nominal pressure used in (MF) equations

%[INERTIA]
MASS    = 9.3           ;% 21 Tire mass
IXX     = 0.391         ;% 22 Tire diametral moment of inertia
IYY     = 0.736         ;% 23 Tire polar moment of inertia
BELT_MASS = 7.247       ;% 24 Belt mass (Swift)
BELT_IXX = 0.3519       ;% 25 Belt diametral moment of inertia (Swift)
BELT_IYY = 0.5698       ;% 26 Belt polar moment of inertia (Swift)
GRAVITY = 0             ;% 27 Gravity acting on belt in Z direction (Swift)

%[VERTICAL]
FNOMIN  = 4000          ;% 28 Nominal wheel load
VERTICAL_STIFFNESS = 209651 ;% 29 Tire vertical stiffness
VERTICAL_DAMPING = 50   ;% 30 Tire vertical damping
MC_CONTOUR_A = 0        ;% 31 Motorcycle contour ellipse A
MC_CONTOUR_B = 0        ;% 32 Motorcycle contour ellipse B
BREFF   = 8.386         ;% 33 Low load stiffness of effective rolling radius
DREFF   = 0.25826       ;% 34 Peak value of effective rolling radius
FREFF   = 0.07394       ;% 35 High load stiffness of effective rolling radius
Q_RE0   = 0.9974        ;% 36 Ratio of Free tire radius with nominal tire radius
Q_V1    = 0.0007742     ;% 37 Tire radius increase with speed
Q_V2    = 0.04667       ;% 38 Vertical stiffness increase with speed
Q_FZ2   = 15.4          ;% 39 Quadratic term in load vs. deflection
Q_FCX   = 0             ;% 40 Longitudinal force influence on vertical stiffness
Q_FCY   = 0             ;% 41 Lateral force influence on vertical stiffness
Q_FCY2  = 0             ;% 42 Explicit load dependency for including the lateral force influence on vertical stiffness
Q_CAM   = 0             ;% 43 Stiffness reduction due to camber
Q_CAM1  = 0             ;% 44 Linear load dependent camber angle influence on vertical stiffness
Q_CAM2  = 0             ;% 45 Quadratic load dependent camber angle influence on vertical stiffness
Q_CAM3  = 0             ;% 46 Linear load and camber angle dependent reduction on vertical stiffness
Q_FYS1  = 0             ;% 47 Combined camber angle and side slip angle effect on vertical stiffness (constant)
Q_FYS2  = 0             ;% 48 Combined camber angle and side slip angle linear effect on vertical stiffness
Q_FYS3  = 0             ;% 49 Combined camber angle and side slip angle quadratic effect on vertical stiffness
Q_FZ3   = 0             ;% 50 Stiffness reduction due to camber (obsolete in 5.2, use Q_CAM instead)
PFZ1    = 0.7098        ;% 51 Pressure effect on vertical stiffness
BOTTOM_OFFST = 0        ;% 52 Distance to rim when bottoming starts to occur
BOTTOM_STIFF = 0        ;% 53 Vertical stiffness of bottomed tire

%[STRUCTURAL]
LONGITUDINAL_STIFFNESS = 358066 ;% 54 Tire overall longitudinal stiffness
LATERAL_STIFFNESS = 102673 ;% 55 Tire overall lateral stiffness
YAW_STIFFNESS = 4793    ;% 56 Tire overall yaw stiffness
FREQ_LONG = 77.17       ;% 57 Undamped frequency fore/aft and vertical mode (Swift)
FREQ_LAT = 42.41        ;% 58 Undamped frequency lateral mode (Swift)
FREQ_YAW = 53.49        ;% 59 Undamped frequency yaw and camber mode (Swift)
FREQ_WINDUP = 58.95     ;% 60 Undamped frequency wind-up mode (Swift)
DAMP_LONG = 0.056       ;% 61 Dimensionless damping fore/aft and vertical mode (Swift)
DAMP_LAT = 0.037        ;% 62 Dimensionless damping lateral mode (Swift)
DAMP_YAW = 0.007        ;% 63 Dimensionless damping yaw and camber mode (Swift)
DAMP_WINDUP = 0.05      ;% 64 Dimensionless damping wind-up mode (Swift)
DAMP_RESIDUAL = 0.002   ;% 65 Residual damping (proportional to stiffness)
DAMP_VLOW = 0.001       ;% 66 Additional low speed damping (proportional to stiffness)
Q_BVX   = 0.364         ;% 67 Load and speed influence on in-plane translation stiffness (Swift)
Q_BVT   = 0.065         ;% 68 Load and speed influence on in-plane rotation stiffness (Swift)
PCFX1   = 0.27504       ;% 69 Tire overall longitudinal stiffness vertical deflection dependency linear term
PCFX2   = 0             ;% 70 Tire overall longitudinal stiffness vertical deflection dependency quadratic term
PCFX3   = 0             ;% 71 Tire overall longitudinal stiffness pressure dependency
PCFY1   = 0.16365       ;% 72 Tire overall lateral stiffness vertical deflection dependency linear term
PCFY2   = 0             ;% 73 Tire overall lateral stiffness vertical deflection dependency quadratic term
PCFY3   = 0.24993       ;% 74 Tire overall lateral stiffness pressure dependency
PCMZ1   = 0             ;% 75 Tire overall yaw stiffness pressure dependency

%[CONTACT_PATCH]
Q_RA1   = 0.671         ;% 76 Square root term in contact length equation (Swift)
Q_RA2   = 0.733         ;% 77 Linear term in contact length equation (Swift)
Q_RB1   = 1.059         ;% 78 Root term in contact width equation (Swift)
Q_RB2   = -1.1878       ;% 79 Linear term in contact width equation (Swift)
ELLIPS_SHIFT = 0.8335   ;% 80 Scaling of distance between front and rear ellipsoid (Swift)
ELLIPS_LENGTH = 1.471   ;% 81 Semimajor axis of ellipsoid (Swift)
ELLIPS_HEIGHT = 0.9622  ;% 82 Semiminor axis of ellipsoid (Swift)
ELLIPS_ORDER = 1.5174   ;% 83 Order of ellipsoid (Swift)
ELLIPS_MAX_STEP = 0.025 ;% 84 Maximum height of road step (Swift)
ELLIPS_NWIDTH = 10      ;% 85 Number of parallel ellipsoids (Swift)
ELLIPS_NLENGTH = 10     ;% 86 Number of ellipsoids at sides of contact patch (Swift)
ENV_C1  = 0             ;% 87 Effective height attenuation
ENV_C2  = 0             ;% 88 Effective plane angle attenuation
Q_A1    = 0             ;% 89 Square root load term in contact length (MF 5.2-6.0 only)
Q_A2    = 0             ;% 90 Linear load term in contact length (MF 5.2-6.0 only)

%[INFLATION_PRESSURE_RANGE]
PRESMIN = 186000        ;% 91 Minimum allowed inflation pressure
PRESMAX = 255000        ;% 92 Maximum allowed inflation pressure

%[VERTICAL_FORCE_RANGE]
FZMIN   = 0             ;% 93 Minimum allowed wheel load
FZMAX   = 10000         ;% 94 Maximum allowed wheel load

%[LONG_SLIP_RANGE]
KPUMIN  = -1            ;% 95 Minimum valid wheel slip
KPUMAX  = 1             ;% 96 Maximum valid wheel slip

%[SLIP_ANGLE_RANGE]
ALPMIN  = -0.96         ;% 97 Minimum valid slip angle
ALPMAX  = 0.96          ;% 98 Maximum valid slip angle

%[INCLINATION_ANGLE_RANGE]
CAMMIN  = -0.105        ;% 99 Minimum valid camber angle
CAMMAX  = 0.105         ;% 100 Maximum valid camber angle

%[SCALING_COEFFICIENTS]
LFZO    = 1             ;% 101 Scale factor of nominal (rated) load
LCX     = 1             ;% 102 Scale factor of Fx shape factor
LMUX    = 1             ;% 103 Scale factor of Fx peak friction coefficient
LEX     = 1             ;% 104 Scale factor of Fx curvature factor
LKX     = 1             ;% 105 Scale factor of slip stiffness
LHX     = 1             ;% 106 Scale factor of Fx horizontal shift
LVX     = 1             ;% 107 Scale factor of Fx vertical shift
LCY     = 1             ;% 108 Scale factor of Fy shape factor
LMUY    = 1             ;% 109 Scale factor of Fy peak friction coefficient
LEY     = 1             ;% 110 Scale factor of Fy curvature factor
LKY     = 1             ;% 111 Scale factor of cornering stiffness
LKYC    = 1             ;% 112 Scale factor of camber stiffness
LKZC    = 1             ;% 113 Scale factor of camber moment stiffness
LHY     = 1             ;% 114 Scale factor of Fy horizontal shift
LVY     = 1             ;% 115 Scale factor of Fy vertical shift
LTR     = 1             ;% 116 Scale factor of peak of pneumatic trail
LRES    = 1             ;% 117 Scale factor for offset of residual torque
LXAL    = 1             ;% 118 Scale factor of alpha influence on Fx
LYKA    = 1             ;% 119 Scale factor of alpha influence on Fx
LVYKA   = 1             ;% 120 Scale factor of kappa induced Fy
LS      = 1             ;% 121 Scale factor of moment arm of Fx
LMX     = 1             ;% 122 Scale factor of overturning moment
LVMX    = 1             ;% 123 Scale factor of Mx vertical shift
LMY     = 1             ;% 124 Scale factor of rolling resistance torque
LMP     = 1             ;% 125 Scale factor of parking moment
LSGKP   = 1             ;% 126 Scale factor of relaxation length of Fx (obsolete after 5.2)
LSGAL   = 1             ;% 127 Scale factor of relaxation length of Fy (obsolete after 5.2)
LMUV    = 0             ;% 128 Friction decay with slip speed (NEW)
LAMU    = 1             ;% 129 Digressive friction factor (NEW)

%[LONGITUDINAL_COEFFICIENTS]
PCX1    = 1.579         ;% 130 Shape factor Cfx for longitudinal force
PDX1    = 1.0422        ;% 131 Longitudinal friction Mux at Fznom
PDX2    = -0.08285      ;% 132 Variation of friction Mux with load
PDX3    = 0             ;% 133 Variation of friction Mux with camber
PEX1    = 0.11113       ;% 134 Longitudinal curvature Efx at Fznom
PEX2    = 0.3143        ;% 135 Variation of curvature Efx with load
PEX3    = 0             ;% 136 Variation of curvature Efx with load squared
PEX4    = 0.001719      ;% 137 Factor in curvature Efx while driving
PKX1    = 21.687        ;% 138 Longitudinal slip stiffness Kfx/Fz at Fznom
PKX2    = 13.728        ;% 139 Variation of slip stiffness Kfx/Fz with load
PKX3    = -0.4098       ;% 140 Exponent in slip stiffness Kfx/Fz with load
PHX1    = 0.00021615    ;% 141 Horizontal shift Shx at Fznom
PHX2    = 0.0011598     ;% 142 Variation of shift Shx with load
PVX1    = 2.0283e-05    ;% 143 Vertical shift Sv/Fz at Fznom
PVX2    = 0.00010568    ;% 144 Variation of shift Sv/Fz with load
RBX1    = 13.046        ;% 145 Slope factor for combined slip Fx reduction
RBX2    = 9.718         ;% 146 Variation of slope Fx reduction with kappa
RBX3    = 0             ;% 147 Influence of camber on stiffness for Fx combined
RCX1    = 0.9995        ;% 148 Shape factor for combined slip Fx reduction
REX1    = -0.4403       ;% 149 Curvature factor of combined Fx
REX2    = -0.4663       ;% 150 Curvature factor of combined Fx with load
RHX1    = -9.968e-05    ;% 151 Shift factor for combined slip Fx reduction
PPX1    = -0.3485       ;% 152 Linear pressure effect on slip stiffness
PPX2    = 0.37824       ;% 153 Quadratic pressure effect on slip stiffness
PPX3    = -0.09603      ;% 154 Linear pressure effect on longitudinal friction
PPX4    = 0.06518       ;% 155 Quadratic pressure effect on longitudinal friction
PTX1    = 0             ;% 156 Relaxation length SigKap0/Fz at Fznom (obsolete after 5.1)
PTX2    = 0             ;% 157 Variation of SigKap0/Fz with load (obsolete after 5.1)
PTX3    = 0             ;% 158 Variation of SigKap0/Fz with exponent of load (obsolete after 5.1)

%[OVERTURNING_COEFFICIENTS]
QSX1    = -0.007764     ;% 159 Overturning moment offset
QSX2    = 1.1915        ;% 160 Camber induced overturning couple
QSX3    = 0.013948      ;% 161 Fy induced overturning couple
QSX4    = 4.912         ;% 162 Mixed load, lateral force and camber on Mx
QSX5    = 1.02          ;% 163 Load effect on Mx with lateral force and camber
QSX6    = 22.83         ;% 164 B-factor of load with Mx
QSX7    = 0.7104        ;% 165 Camber with load on Mx
QSX8    = -0.023393     ;% 166 Lateral force with load on Mx
QSX9    = 0.6581        ;% 167 B-factor of lateral force with load on Mx
QSX10   = 0.2824        ;% 168 Vertical force with camber on Mx
QSX11   = 5.349         ;% 169 B-factor of vertical force with camber on Mx
QSX12   = 0             ;% 170 Camber squared induced overturning moment
QSX13   = 0             ;% 171 Lateral force induced overturning moment
QSX14   = 0             ;% 172 Lateral force induced overturning moment with camber
PPMX1   = 0             ;% 173 Influence of inflation pressure on overturning moment

%[LATERAL_COEFFICIENTS]
PCY1    = 1.338         ;% 174 Shape factor Cfy for lateral forces
PDY1    = 0.8785        ;% 175 Lateral friction Muy
PDY2    = -0.06452      ;% 176 Variation of friction Muy with load
PDY3    = 0             ;% 177 Variation of friction Muy with squared camber
PEY1    = -0.8057       ;% 178 Lateral curvature Efy at Fznom
PEY2    = -0.6046       ;% 179 Variation of curvature Efy with load
PEY3    = 0.09854       ;% 180 Zero order camber dependency of curvature Efy
PEY4    = -6.697        ;% 181 Variation of curvature Efy with camber
PEY5    = 0             ;% 182 Camber curvature Efc
PKY1    = -15.324       ;% 183 Maximum value of stiffness Kfy/Fznom
PKY2    = 1.715         ;% 184 Load at which Kfy reaches maximum value
PKY3    = 0.3695        ;% 185 Variation of Kfy/Fznom with camber
PKY4    = 2.0005        ;% 186 Curvature of stiffness Kfy
PKY5    = 0             ;% 187 Peak stiffness variation with camber squared
PKY6    = -0.8987       ;% 188 Camber stiffness factor
PKY7    = -0.23303      ;% 189 Load dependency of camber stiffness factor
PHY1    = -0.001806     ;% 190 Horizontal shift Shy at Fznom
PHY2    = 0.00352       ;% 191 Variation of shift Shy with load
PHY3    = 0             ;% 192 Variation of shift Shy with camber (MF 5.2 only)
PVY1    = -0.00661      ;% 193 Vertical shift in Svy/Fz at Fznom
PVY2    = 0.03592       ;% 194 Variation of shift Svy/Fz with load
PVY3    = -0.162        ;% 195 Variation of shift Svy/Fz with camber
PVY4    = -0.4864       ;% 196 Variation of shift Svy/Fz with camber and load
RBY1    = 10.622        ;% 197 Slope factor for combined Fy reduction
RBY2    = 7.82          ;% 198 Variation of slope Fy reduction with alpha
RBY3    = 0.002037      ;% 199 Shift term for alpha in slope Fy reduction
RBY4    = 0             ;% 200 Influence of camber on stiffness of Fy combined
RCY1    = 1.0587        ;% 201 Shape factor for combined Fy reduction
REY1    = 0.3148        ;% 202 Curvature factor of combined Fy
REY2    = 0.004867      ;% 203 Curvature factor of combined Fy with load
RHY1    = 0.009472      ;% 204 Shift factor for combined Fy reduction
RHY2    = 0.009754      ;% 205 Shift factor for combined Fy reduction with load
RVY1    = 0.05187       ;% 206 Kappa induced side force Svyk/Muy*Fz at Fznom
RVY2    = 0.0004853     ;% 207 Variation of Svyk/Muy*Fz with load
RVY3    = 0             ;% 208 Variation of Svyk/Muy*Fz with camber
RVY4    = 94.63         ;% 209 Variation of Svyk/Muy*Fz with alpha
RVY5    = 1.8914        ;% 210 Variation of Svyk/Muy*Fz with kappa
RVY6    = 23.8          ;% 211 Variation of Svyk/Muy*Fz with atan(kappa)
PPY1    = -0.6255       ;% 212 Pressure effect on cornering stiffness magnitude
PPY2    = -0.06523      ;% 213 Pressure effect on location of cornering stiffness peak
PPY3    = -0.16666      ;% 214 Linear pressure effect on lateral friction
PPY4    = 0.2811        ;% 215 Quadratic pressure effect on lateral friction
PPY5    = 0             ;% 216 Influence of inflation pressure on camber stiffness
PTY1    = 0             ;% 217 Peak value of relaxation length SigAlp0/R0
PTY2    = 0             ;% 218 Value of Fz/Fznom where SigAlp0 is extreme

%[ROLLING_COEFFICIENTS]
QSY1    = 0.00702       ;% 219 Rolling resistance torque coefficient
QSY2    = 0             ;% 220 Rolling resistance torque depending on Fx
QSY3    = 0.001515      ;% 221 Rolling resistance torque depending on speed
QSY4    = 8.514e-05     ;% 222 Rolling resistance torque depending on speed^4
QSY5    = 0             ;% 223 Rolling resistance torque depending on camber squared
QSY6    = 0             ;% 224 Rolling resistance torque depending on load and camber squared
QSY7    = 0.9008        ;% 225 Rolling resistance torque coefficient load dependency
QSY8    = -0.4089       ;% 226 Rolling resistance torque coefficient pressure dependency

%[ALIGNING_COEFFICIENTS]
QBZ1    = 12.035        ;% 227 Trail slope factor for trail Bpt at Fznom
QBZ2    = -1.33         ;% 228 Variation of slope Bpt with load
QBZ3    = 0             ;% 229 Variation of slope Bpt with load squared
QBZ4    = 0.176         ;% 230 Variation of slope Bpt with camber
QBZ5    = -0.14853      ;% 231 Variation of slope Bpt with absolute camber
QBZ9    = 34.5          ;% 232 Slope factor Br of residual torque Mzr
QBZ10   = 0             ;% 233 Slope factor Br of residual torque Mzr
QCZ1    = 1.2923        ;% 234 Shape factor Cpt for pneumatic trail
QDZ1    = 0.09068       ;% 235 Peak trail DptPrime = Dpt*(Fz/Fznom*R0)
QDZ2    = -0.00565      ;% 236 Variation of peak Dpt with load
QDZ3    = 0.3778        ;% 237 Variation of peak Dpt with camber
QDZ4    = 0             ;% 238 Variation of peak Dpt with camber squared
QDZ6    = 0.0017015     ;% 239 Peak residual torque Dmr = Dmr/(Fz*R0)
QDZ7    = -0.002091     ;% 240 Variation of peak factor Dmr with load
QDZ8    = -0.1428       ;% 241 Variation of peak factor Dmr with camber
QDZ9    = 0.00915       ;% 242 Variation of peak factor Dmr with camber and load
QDZ10   = 0             ;% 243 Variation of peak factor Dmr with camber squared
QDZ11   = 0             ;% 244 Variation of Dmr with camber squared and load
QEZ1    = -1.7924       ;% 245 Trail curvature Ept at Fznom
QEZ2    = 0.8975        ;% 246 Variation of curvature Ept with load
QEZ3    = 0             ;% 247 Variation of curvature Ept with load squared
QEZ4    = 0.2895        ;% 248 Variation of curvature Ept with sign of Alpha-t
QEZ5    = -0.6786       ;% 249 Variation of Ept with camber and sign Alpha-t
QHZ1    = 0.0014333     ;% 250 Trail horizontal shift Sht at Fznom
QHZ2    = 0.0024087     ;% 251 Variation of shift Sht with load
QHZ3    = 0.24973       ;% 252 Variation of shift Sht with camber
QHZ4    = -0.21205      ;% 253 Variation of shift Sht with camber and load
SSZ1    = 0.00918       ;% 254 Nominal value of s/R0: effect of Fx on Mz
SSZ2    = 0.03869       ;% 255 Variation of distance s/R0 with Fy/Fznom
SSZ3    = 0             ;% 256 Variation of distance s/R0 with camber
SSZ4    = 0             ;% 257 Variation of distance s/R0 with load and camber
PPZ1    = -0.4408       ;% 258 Linear pressure effect on pneumatic trail
PPZ2    = 0             ;% 259 Influence of inflation pressure on residual aligning torque

%[TURNSLIP_COEFFICIENTS]
PDXP1   = 0.4           ;% 260 Peak Fx reduction due to spin parameter
PDXP2   = 0             ;% 261 Peak Fx reduction due to spin with varying load parameter
PDXP3   = 0             ;% 262 Peak Fx reduction due to spin with kappa parameter
PKYP1   = 1             ;% 263 Cornering stiffness reduction due to spin
PDYP1   = 0.4           ;% 264 Peak Fy reduction due to spin parameter
PDYP2   = 0             ;% 265 Peak Fy reduction due to spin with varying load parameter
PDYP3   = 0             ;% 266 Peak Fy reduction due to spin with alpha parameter
PDYP4   = 0             ;% 267 Peak Fy reduction due to square root of spin parameter
PHYP1   = 1             ;% 268 Fy-alpha curve lateral shift limitation
PHYP2   = 0.15          ;% 269 Fy-alpha curve maximum lateral shift parameter
PHYP3   = 0             ;% 270 Fy-alpha curve maximum lateral shift varying with load parameter
PHYP4   = -4            ;% 271 Fy-alpha curve maximum lateral shift parameter
PECP1   = 0.5           ;% 272 Camber w.r.t. spin reduction factor parameter in camber stiffness
PECP2   = 0             ;% 273 Camber w.r.t. spin reduction factor varying with load parameter in camber stiffness
QDTP1   = 10            ;% 274 Pneumatic trail reduction factor due to turn slip parameter
QCRP1   = 0.2           ;% 275 Turning moment at constant turning and zero forward speed parameter
QCRP2   = 0.1           ;% 276 Turn slip moment (at alpha=90deg) parameter for increase with spin
QBRP1   = 0.1           ;% 277 Residual (spin) torque reduction factor parameter due to side slip
QDRP1   = 1             ;% 278 Turn slip moment peak magnitude parameter
QDRP2   = 0             ;% 279 Turn slip moment peak position parameter

end

properties(Access = private)
  
% Computed parameters that depend on other parameters;
% computed by on_params_changed()
Q_A1_   = 0;
Q_A2_   = 0;
Q_FZ1_  = 0;

end

properties(Constant)

% get persistent copy of param_info table
PARAM_INFO = init_param_info();
end

end % end classdef

%% [Eqn (4.E6a) Page 178 - Book]
% Apply an epsilon so that returned x is
% slightly larger than x in magnitude and never zero.
function x = make_nonzero(x)
  epsilon = 1e-6;
  if x < 0
    x = x - epsilon;
  else
    x = x + epsilon;
  end
  %x = x + ((x >= 0) * 2 - 1)*epsilon;
end

%% helper functions

function n = vec3_norm(v)
  n = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
end

function v = vec3_normalize(v)
  n = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
  % norm assumed to be > 0
  %assert(n > 0);
  v = v / n;
end

function d = vec3_dot(a, b)
  d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
end

function c = vec3_cross(a, b)
  c = [a(2)*b(3) - a(3)*b(2);
       a(3)*b(1) - a(1)*b(3);
       a(1)*b(2) - a(2)*b(1)];
end

% Fortran-style SGN function: returns -1 if (x < 0) and 1 if (x >=) 0.
function s = SGN(x)
  %s = ((x >= 0) * 2 - 1);  % vectorized version
  s = 1;
  if (x < 0) %then
    s = -1;
  end
end

% Square a value
function y = SQR(x)
  y = x*x;
end

% Raise x to the power of e
function y = POW(x, e)
  y = x ^ e;
end

% clamp to value between min_val and max_val
function clamped = CLAMP(x, min_val, max_val)
  clamped = min(max(x, min_val), max_val);
end

% smoothstep function for smooth transition between s = 0 and s = 1
% A cubic with first derivatives == 0 at ends, y = 0 at 0, and y = 1 at 1.
% smootherstep is similar, but both first and second derivatives are 0 at ends.
% smoothstep is closest to 1-2*cos(pi*s); smootherstep is more gradual at ends,
% steeper at s = 0.5.  Both are much more efficient than 1-2*cos(pi*s).
% Equivalent of sin(s*pi/2) is 2*SMOOTHSTEP((s+1)/2)-1
function y = SMOOTHSTEP(s)
  y = s*s*(3 - 2*s);
end

% Similar to sin(s*pi/2)
% Derivative at s==0 is 1, derivative at s==1 is 0.
function y = SHARPSTEP(s)
  y = 2*SMOOTHSTEP((s + 1)/2) - 1;
end

function y = SMOOTHERSTEP(s)
  y = s*s*s*(s*(s * 6 - 15) + 10);
end

% Linear interpolation
function y = LERP(x, x0, x1, y0, y1)
  y = y0 + (x - x0) / (x1 - x0) * (y1 - y0);
end

%% Sine form of magic formula
% y = D*sin(C*atan(B*x - E*(B*x - atan(B*x))));
function y = sin_magic_formula(x, B, C, D, E)
    a = B*x;
    %assert(E <= 1);
    if (E ~= 0) %then
      a = a - E*(a - atan(a));
    end
    y = D * sin(C * atan(a));
end

%% Sine form of magic formula with E == 0
function y = sin_magic_formula_3(x, B, C, D)
    y = D * sin(C * atan(B*x));
end

%% Cosine form of magic formula
% y = D*cos(C*atan(B*x - E*(B*x - atan(B*x))));
function y = cos_magic_formula(x, B, C, D, E)
    a = B*x;
    if (E ~= 0) %then
    %assert(E <= 1);
      a = a - E*(a - atan(a));
    end
    y = D * cos(C * atan(a));
end

%% Cosine form of magic formula with E == 0
function y = cos_magic_formula_3(x, B, C, D)
    y = D * cos(C * atan(B*x));
end

%% Cosine form of magic formula with C,D = 1 and E == 0
function y = cos_magic_formula_1(x, B)
    y = cos(atan(B*x));
end

%% Parallel trig and magic formula functions
% (two for the price of one in C++)
%
function [y1,y2] = two_sin(x1, x2)
  y1 = sin(x1);
  y2 = sin(x2);
end

function [y1,y2] = two_cos(x1, x2)
  y1 = cos(x1);
  y2 = cos(x2);
end

function [y1,y2] = two_tan(x1, x2)
  y1 = tan(x1);
  y2 = tan(x2);
end

function [x1,x2] = two_atan(y1, y2)
  x1 = atan(y1);
  x2 = atan(y2);
end

function [y1,y2] = two_sin_magic_formula(x1, x2, B1, B2, C1, C2, D1, D2, E1, E2)
    a1 = B1*x1;
    a2 = B2*x2;
    if (E1 ~= 0 || E2 ~= 0) %then
      %assert(E1 <= 1 && E2 <= 0);
      [a_a1,a_a2] = two_atan(a1,a2);
      a1 = a1 - E1*(a1 - a_a1);
      a2 = a2 - E2*(a2 - a_a2);
    end
    [a_a1,a_a2] = two_atan(a1,a2);
    [y1,y2] = two_sin(C1*a_a1, C2*a_a2);
    y1 = D1 * y1;
    y2 = D2 * y2;
end

function [y1,y2] = two_cos_magic_formula(x1, x2, B1, B2, C1, C2, D1, D2, E1, E2)
    a1 = B1*x1;
    a2 = B2*x2;
    if (E1 ~= 0 || E2 ~= 0) %then
      %assert(E1 <= 1 && E2 <= 0);
      [a_a1,a_a2] = two_atan(a1,a2);
      a1 = a1 - E1*(a1 - a_a1);
      a2 = a2 - E2*(a2 - a_a2);
    end
    if (E2 ~= 0) %then
    end
    [a_a1,a_a2] = two_atan(a1,a2);
    [y1,y2] = two_cos(C1*a_a1, C2*a_a2);
    y1 = D1 * y1;
    y2 = D2 * y2;
end

%% A simple one-dimensional solver using the secant method.
%
% Finds x where |f(x)| < ytol.
% x0 and x1 are two starting guesses,
% which ideally are above and below the solution.
% Returns x = [] if no solution could be found.

% Makes a minimum of 3 calls to the function to be zeroed;
% typically 4 or 5.
function x = solve_secant(fun, x0, x1, ytol)

  % Assume that the first two guesses are wrong
  y0 = fun(x0);
  y1 = fun(x1);

  max_calls = 8;
  calls = 2;
  while (true) %then
    
    % interpolate (or extrapolate) to y == 0
    x = x1 - (x1 - x0) / (y1 - y0) * y1;
    
    % if (x1 - x0) is zero then fun() is not continuous,
    % possibly due to roundoff issues in fun().
    % if (y1 - y0) is zero, then we have arrived
    % at a local minimum of fun().
    %assert(isfinite(x) && x ~= x1);
    
    % evaluate new guess and check convergence
    y = fun(x);
    
    calls = calls + 1;
    if (max(abs(y)) <= ytol) %then
      return;
    end
%     if (calls >= 8) %then
%       fprintf('solve_secant: slow convergence: %d calls, y = %.2f\n', calls, abs(y)/ytol);
%     end
    if (calls == max_calls) %then
      %assert(false); % increase max_calls or choose better initial guess
      break;
    end
    
    % Keep the best in x0
    if (abs(y1) < abs(y0)) %then
      x0 = x1;
      y0 = y1;
    end
    
    % replace x1,y1 for next iteration
    x1 = x;
    y1 = y;
  end
  % can't converge: try increasing max_calls
  %assert(false);
end

function J = one_sided_jacobian(f, x, y)
  m = length(y);
  n = length(x);
  J = zeros(m, n);
  delta = sqrt(eps);
  for j=1:n
    xj = x(j);
    h = max(delta,delta*abs(xj));
    if h == 0, h = delta; end
    x(j) = xj + h;
    J(:,j) = (f(x) - y) / (x(j) - xj);
    x(j) = xj;
  end
end
