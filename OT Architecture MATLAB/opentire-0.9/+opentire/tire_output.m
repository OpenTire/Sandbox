classdef tire_output
%TIRE_OUTPUT Tire output base class

properties
right_tire = 0;           % [bool]

force = [nan;nan;nan];    % Fx, Fy, Fz
moment = [nan;nan;nan];   % Mx, My, Mz

tire_contact  =     [nan;nan;nan];
tire_contact_vel =  [nan;nan;nan];
tire_slip_vel =     [nan;nan];

wheel_omega     = nan;    % [rad/s] wheel angular velocity
rolling_radius  = nan;    % [m] effective rolling radius

kappa = nan;              % [-]
alpha = nan;              % [rad]
gamma = nan;              % [rad]
path_curvature = nan;     % [1/m] 1/R
phi = nan;                % [1/m] effective turnslip (path curvature plus gamma-induced)

loaded_radius = nan;      % [m]
rho = nan;                % [m]
vertical_stiffness = nan; % [N/m]
vertical_damping = nan;   % [N/(m/s)]

trail = [nan;nan];        % [m] longitudinal and lateral pneumatic trails
grip  = [nan;nan];        % [-] lon and lat grip

relaxation     = [nan;nan];         % [m] lon and lat relaxation lengths
slip_stiffness = [nan;nan;nan;nan]; % nominal slip stiffnesses (dFx/dkappa, dFy/dalpha, dFy/dgamma, dMz/dphi)
contact_size   = [nan;nan];         % [m] X and Y contact patch dimensions (full size, not half)

W = nan(3,3);             % Rotation to inertial frame from Tydex W frame
R = nan(3,3);             % Rotation to inertial frame from wheel frame

end % properties

methods
  
function output = change_coordinate_system(output, from, to)
  if (from > to)
    t = from; from = to; to = t;
  end
  if (from == 0)
    if (to == 2)
      output.force(2) = -output.force(2);
      output.moment(3) = -output.moment(3);
      output.phi = -output.phi;
      output.path_curvature = -output.path_curvature;
      output.tire_slip_vel(2) = -output.tire_slip_vel(2);
      output.tire_contact_vel(2) = -output.tire_contact_vel(2);
      output.W(:,2) = -output.W(:,2);
      output.R(:,2) = -output.R(:,2);
    end
  else
    assert(false); % not yet supported
  end
end
  
end % methods

end % classdef
