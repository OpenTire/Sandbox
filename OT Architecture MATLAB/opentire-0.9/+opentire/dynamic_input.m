classdef dynamic_input
%DYNAMIC_INPUT Base class for compute_dynamic() input data
properties
  
right_tire = false;           % [bool] 1 if mounted on right, 0 if left

time = 0;                     % [s] time since simulation start (optional)
distance = 0;                 % [m] simulated distance traveled (optional)

frame_omega = [nan;nan;nan];  % [rad/s] angular velocity of rotating frame
wheel_center = [nan;nan;nan]; % [m] wheel center in inertial frame
wheel_center_vel = [nan;nan;nan]; % [m/s] wheel center velocity in inertial frame
wheel_axis = [nan;nan;nan];   % wheel axis normal in inertial frame
wheel_omega = nan;            % [rad/s] wheel angular velocity about wheel_axis

grip_adjust = [1;1];          % [-] grip adjustment
pressure = 0;                 % [pascal] inflation pressure

end % properties
end % classdef

