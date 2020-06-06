%% opentire.brush_tire
%
% This class implements the brush model described in Chapter 3 of
% "Tire and Vehicle Dynamics, 3rd Edition", by Hans Pacejka.
%
% Based on TreadSim.m, written by H.B. Pacejka.
%
classdef (Sealed = true) brush_tire < opentire.tire

methods(Access = public)

%% Constructor
% Construct from tire filename or existing handle
function tire = brush_tire(params)
  if nargin == 0
    params = 0;
  end
  
  if (ischar(params)) %then
    
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
    warning('brush_tire: Unsupported constructor parameter');
    tire = []; % return empty handle to indicate error
  end
end

%% Create options structure
function options = new_options(tire)
  options.log_deflections = false;
end

%% Create input structure with nominal values for steady-state simulations
function input = new_steady_input(tire, right_tire)
  
  if (nargin < 2)  %then
    right_tire = false;
  end
  
  input = opentire.steady_input();
  
  input.right_tire = right_tire;
  input.Fz = tire.Fz;
  input.Vx = 100/3.6;
  input.loaded_radius = nan;
  input.kappa = 0;
  input.alpha = 0;
  input.gamma = 0;
  input.path_curvature = 0;
  input.pressure = 0;%tire.INFLPRES;
end

function input = new_dynamic_input(tire, right_tire)
  if (nargin < 2)  %then
    right_tire = false;
  end
  
  input = opentire.dynamic_input();
  
  input.right_tire = right_tire;
end

%% compute_steady, compute_dynamic -- Tire model computation functions
function [output, states_dot, algebraic_loops] = compute_dynamic(t, input, road, states, algebraic_loops, options)
  states_dot = [];
  output = opentire.brush_tire_output;
  assert(false); % not yet implemented
end

function output = compute_steady(tire, input, options)
  output = opentire.brush_tire_output;
  
  Fz = input.Fz;
  Vx = input.Vx;
  kappa = input.kappa;
  alpha = input.alpha;
  gamma = input.gamma;
  phi_t = input.path_curvature;
  
  output.kappa = kappa;
  output.alpha = alpha;
  output.gamma = gamma;
  output.path_curvature = phi_t;
  
  output.rho = Fz/tire.vertical_stiffness;
  output.loaded_radius = tire.unloaded_radius - output.rho;
  output.vertical_stiffness = tire.vertical_stiffness;
  output.vertical_damping = 0;
  
  output.rolling_radius = tire.Re;
  
  Vcx = Vx;
  Vsx = Vcx * input.kappa;
  Vsy = Vcx * tan(input.alpha);
  Vr    = Vcx - Vsx;
  output.wheel_omega = Vr / tire.Re;
  
  output.contact_size = tire.contact_size;
  
  % REVIEW: This common code can be output method
  sin_alpha = sin(alpha); cos_alpha = cos(alpha);
  sin_gamma = sin(gamma); cos_gamma = cos(gamma);
  
  wheel_axis = [sin_alpha*cos_gamma; cos_alpha*cos_gamma; sin_gamma];
  road_normal = [0;0;1];
  
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
    
  output.R = R;
  output.W = W;
  
  [Fx,Fy,Mz,deflections] = tire.brush_model(Fz, Vx, kappa, alpha, gamma, phi_t, options.log_deflections);
  
  if (~isempty(deflections))
    output.p = deflections.p;
    output.e = deflections.e;
    output.b = deflections.b;
  end
  
  output.force = [Fx,Fy,Fz];
  output.moment = [0;0;Mz];
  
end

function [Fx,Fy,Mz,out] = brush_model(tire, Fz, Vx, kappa, alpha, gamma, curvature, log_deflections)
  
  rows = tire.rows; %#ok<*PROPLC>
  n = tire.elements;
  
  out = [];
  if log_deflections
    out.p = zeros(3,rows,n+1);
    out.e = zeros(2,rows,n+1);
    out.b = zeros(2,rows,n+1);
  end
  
  a = tire.contact_size(1)/2;
  b = tire.contact_size(2)/2;
  
  rolling_radius = tire.Re;
  
  CFkappa = tire.CFkappaFz * Fz; %[N/-]
  % from this: cp = (1/(nrow*2*a*a)) * CFkappa;
  
  % resulting stiffness of tread elements,
  % per row, per unit of circumference
  Cp = CFkappa * (1 / (rows*2*a*a));
  %theta = (1/3) * CFkappa / (mu0*Fz); % not used
  
  gamma = gamma + tire.conicity_gamma;
  sin_gamma = sin(gamma);
  
  cs0 = tire.ply_steer_alpha;                       % initial belt slope due to ply-steer
  cc0 = (1/rolling_radius) * (1-tire.eps_gamma) * sin_gamma; % initial belt lateral curvature due to camber and conicity
  
  % initial belt lateral position at contact center
  % resulting from camber at mu = 0. This is a guess.
  % Then, y0 = cs0*a at leading edge (x = a)
  ym0 = -cc0*a*a/2;
  
  % moment arm of Fx due to sideways rolling at camber
  y_gamma_roll = tire.eps_y_gamma_roll * sin_gamma * b;
  if abs(y_gamma_roll) > b
    y_gamma_roll = b * sign(gamma);
  end
  
  y_rows = 0;
  if rows > 1
    y_rows = linspace(-0.5,0.5,rows) * tire.contact_size(2) * tire.effective_width_scale;
  end
  
  psi_dot = -Vx*curvature; % curvature = 1/R
  Vcx = Vx * cos(alpha);
  Vsx = -Vcx * kappa;
  Vsy = -Vcx * tan(alpha);
  
  % zero slip velocities are a problem: just make them small
  if (Vsx == 0)
    Vsx = 1e-5;
  end

  Vr    = Vcx - Vsx;
  omega = Vr / rolling_radius;
  
  Fy = 0; % algebraic loop variables
  Mz_prime = 0;
  max_carcass_iter = 35;
  for carcass_iter = 1:max_carcass_iter
    
    Fy_alg_loop = Fy;
    Mz_alg_loop = Mz_prime;
    
    % compute coefficients of a quadratic in bx
    % that yield carcass deflection at bx
    %
    % by = y_row + ym + bx*cs + cc*bx*bx/2;
    
    cs = cs0;
    cc = cc0;
    ym = ym0;
    y_Fyroll = 0;
    if (~tire.rigid)
      
      slope = Mz_prime / tire.C_yaw;
      cs    = cs0 + slope;         % belt slope at contact center w.r.t. wheel plane

      curve = Fy / tire.C_bend;
      cc    = cc0 - curve;         % belt lateral curvature in contact zone

      ydefl = Fy / tire.C_lat;
      ym    = ym0 + ydefl;         % belt lateral position at contact center w.r.t. wheel plane
    
      y_Fyroll = -tire.eps_y_Fy * ydefl;
    end
    
    % effective Mz moment arm for Fx
    ym_eff = ym + y_gamma_roll + y_Fyroll;

    % integration sums
    integration_sums = [0;0;0];
    
    for irow = 1:rows
      y_row  = y_rows(irow);
      
      % velocities for this row
      %Vcx_row = Vcx - y_row * psi_dot;
      
      Vsx_row = Vsx - y_row * (psi_dot + omega * (1-tire.eps_Re_gamma) * sin_gamma);
      Vsy_row = Vsy;
      Vr_row  = Vr;
      
      delta_x = 2*a / n;
      delta_t = delta_x / (Vr_row + 1e-3*SGN(Vr_row));  % avoid division by zero
      
      xb = a;  % from a to -a, step -delta_x
      e = [0;0];
      p = [0;0;0];
      
      % store first point
      if log_deflections
        out.p(:,irow,1) = p;
        out.e(:,irow,1) = e;
        yb = y_row + ym + cs*xb + cc*xb*xb/2;
        out.b(:,irow,1) = [xb;yb];
      end
      
      % **  begin i-loop: passage through contact length
      
      for i = 2:n+1
        % (xb,yb) = point B, including carcass deflection
        
        xb = xb - delta_x; % or xb = a - (i-1)*delta_x
        
        yb = y_row + ym + cs*xb + cc*xb*xb/2;
        
        % effective moment arm
        yb_eff = (ym_eff - ym) + yb;

        % Velocities of belt at the midpoint between
        % point B and point B of previous timestep
        xb_mid = (xb + delta_x/2);
        dyb_dxb_mid = cs + cc*xb_mid;
        Vb = [Vsx_row
              Vsy_row + psi_dot*xb_mid - Vr_row*dyb_dxb_mid];
        
        norm_Vb = norm(Vb);
        
        % Speed-adjusted grip
        mu = tire.mu0 / (1 + tire.grip_speed_adjust * norm_Vb);
        
        % Vertical pressure distribution
        p(3) = Fz * 3*(a*a - xb*xb) / (4*a*a * a*rows);
        
        delta_s = delta_t * Vb;
        
        if i == 2
          % first time through, use tiny displacement
          % in opposite direction of velocity vector
          tiny_e = -1e-5;
          e = tiny_e * (Vb ./ norm_Vb);
        end

        % grip-limited pressure limit
        p_limit = mu * p(3);
        
        % grip-limited displacement limit
        e_limit = p_limit / Cp;
          
        % compute new displacement assuming element
        % is sticking to road (no sliding)
        e_ns = e - delta_s;
        
        if i == 2
          % see if we should start out sliding or not
          sliding = (dot(e_ns,e_ns) > e_limit^2);
        end

        if ~sliding
          
          % update e
          e = e_ns;

          sliding = (dot(e,e) >= (e_limit - 1e-6)^2);
          if sliding
            % clamp to avoid overshoot
            e = e .* (e_limit / norm(e));
          end
          
        else % sliding

          %norm_e = norm(e);
          %g = norm_e * (dot(e_ns,e_ns) - dot(e,e)) / (2 * dot(e,e_ns));
          %e = (1 - g/norm_e) .* e - delta_s;

          % g * norm_e is sliding distance
          g = (dot(e_ns,e_ns) - dot(e,e)) / (2 * dot(e,e_ns));
          
          % update e
          e = (1 - g) .* e - delta_s;
          
          % at low velocities with large time steps,
          % the computation above can over- or undershoot e_limit
          % If we were sliding, then set magnitude of e to e_limit.
          e = e .* (e_limit / norm(e));
          
          sliding = (g >= 0);
        end
        
        p(1:2) = e * Cp;
        
        % Integrate row forces and moments
        integration_sums = integration_sums + [p(1:2); xb*p(2) - yb_eff*p(1)];
        
        if log_deflections
          out.p(:,irow,i) = p;
          out.e(:,irow,i) = e;
          out.b(:,irow,i) = [xb;yb];
        end
      end  % **  end i-loop (passage through contact length)
    end  % **  end i_row-loop (one, two or three rows of elements)
    
    integration_sums = integration_sums .* 2*a/n;
    Fx = integration_sums(1);
    Fy = integration_sums(2);
    Mz = integration_sums(3);
    
    yc = 0;
    if tire.rigid
      yc = Fy * tire.C;
      % If rigid, adjust y coordinates by
      % moment arm adjustment
      if log_deflections
        out.b(2,:,:) = out.b(2,:,:) + yc;
      end
    end
    
    Mz = Mz - Fx * (tire.y0 - yc); % (yc == 0 if rigid)
    
    Mz_prime = Mz + Fx * ym_eff * tire.eps_y_prime;      % for torsion (slope) calculation
    
    % Exit when algebraic loop has been closed (i.e. inputs same as outputs)
    Fy_tol = tire.Fz * 1e-2;
    Mz_tol = 50 * 1e-2;
    if (tire.rigid || abs(Fy - Fy_alg_loop) < Fy_tol && abs(Mz_prime - Mz_alg_loop) < Mz_tol)
      break;
    end
  end  % carcass deflection iteration
  
  %assert(carcass_iter ~= max_carcass_iter);
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
  tire.model_version = 100;
  tire.right_tire = false;
end
end % methods(Access = protected)

properties
% REVIEW: fix these names...

  rigid  = false;  % 1 if carcass is rigid, else 0
  rows   = 3;     % number of rows (1 or 2 or 3)
  elements = 30;  % number of bristle elements in longitudinal direction
  
  % loaded radius model parameters
  vertical_stiffness = 25*9.81/(1/1000);
  unloaded_radius = 0.300;
  
  %
  Fz     = 3000;   % Nominal normal load [N]
  Re     = 0.300;  % effective rolling radius [m]
  
  contact_size = 2*[0.100 0.080];
  
  effective_width_scale = 0.625;  % effective contact width (scale of contact_size(2))
  
  mu0    = 1.0;    % friction coefficient [-]
  
  grip_speed_adjust = 0.03;   % speed dependency coefficient for mu [s/m] ( mu = mu0/(1+grip_speed_adjust*norm_Vb) );
  
  CFkappaFz = 15;  % long. slip stiffness over Fz [-]
  
  ply_steer_alpha = 0*pi/180;  % ply-steer equivalent slip angle [rad]
  conicity_gamma  = 0*pi/180;  % conicity equivalent camber angle [rad]
  
  % Carcass stiffnesses not used if rigid
  C_yaw  = 6000; % [Nm/rad]
  C_bend = 4000; % [mN]
  C_lat  = 3000/(0.15*0.100); % Fz/(0.15*contact_size(1)/2);  % [N/m]

  % only used if rigid is true
  C      = 0.25/1e5; % c = (1-eps_y_Fy)/C_lat, Mz = Mz - y0*Fx - c*Fx*Fy
  
  %**  correction parameters
  eps_y_prime = 0.5;  % reduction factor for moment (Mz_prime) arm of Fx causing yaw distortion (slope);
  % arm  = (1-eps_y_prime)*ym_eff   (0.5)
  % where ym_eff represents effective lateral displacement of belt in contact zone
  
  eps_y_gamma_roll = 4;   % moment (Mz) arm for Fx due to sideways rolling caused by camber
  % y_gamma_roll = eps_y_gamma_roll * gamma*by; if abs(y_gamma_roll)>b: y_gamma_roll = b*sign(gamma)
  
  eps_y_Fy = 0.0;     % reduction factor moment (Mz) arm for Fx due to sideways rolling caused by Fy
  % and the counter effect of longitudinal deflection
  % total moment {Mz) arm =  ym_eff =  ym0 + ydefl(1-eps_y_Fy) + y_gamma_roll
  
  eps_Re_gamma = 0.0;   % coefficient for reduced change eff. rolling radius left/right caused by camber
  % delta_re_Right == -(1-eps_Re_gamma)*gamma*y_row; suggestion: eps_Re_gamma = eps_xgamma = eps_gamma
  
  eps_gamma = 0.0;   % reduced camber curvature coefficient (will be >0) due to distortion (=eps_ygamma)
  % resulting camber contactline curvature (at mu = 0)
  % 1/Rgamma = (1/re)*(1-eps_gamma)*sin(gamma)
  
  y0 = 0.00;        % initial lateral offset (v0), if rigid: Mz = Mz - c*Fx*Fy - y0*Fx  (0.005)
end

properties(Constant)
% get persistent copy of param_info table
PARAM_INFO = [];
end

end % end classdef

% Fortran-style SGN function: returns -1 if (x < 0) and 1 if (x >=) 0.
function s = SGN(x)
  %s = ((x >= 0) * 2 - 1);  % vectorized version
  s = 1;
  if (x < 0) %then
    s = -1;
  end
end

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
