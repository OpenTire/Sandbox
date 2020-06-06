classdef steady_input
%STEADY_INPUT Base class for compute_steady() input data
properties
  
  right_tire  = false;  % [bool] 1 if tire mounted on right, 0 left.
  
  Fz          = nan;    % [N] Vertical force in N (0 if not known)
  loaded_radius = nan;  % [m] loaded_radius (nan/0 if not known)
  path_curvature = nan; % [1/m] 1/R
  Vx          = nan;    % [m/s]
  wheel_omega = nan;    % [rad/s] (nan/0 if not known)
  kappa       = nan;    % [-] (nan/0 if not known)
  alpha       = nan;    % [rad]
  gamma       = nan;    % [rad]
  pressure    = 0;      % [pascal]  (0 for nominal/unknown)

end % properties

methods
  
function input = change_coordinate_system(input, from, to)
  assert(false); % not yet implemented
end
  
end

end % classdef

