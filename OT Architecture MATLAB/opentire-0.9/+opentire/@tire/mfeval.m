%% diff - Compare two tire objects
% Prints differing parameters to console window.
%
% part of mftire 1.1.0
%function result = diff(tire, compare_tire)
function OUTPUT = mfeval(tire, INPUT, options)
  if (nargin < 3 || isempty(options)) %then
    options = tire.new_options();
  end
  
  [rows, cols] = size(INPUT);
  
  input = tire.new_steady_input;
  
  input.right_tire = false; % or maybe tire.is_right_tire()?
  
  % default values in case not present
  input.pressure = 0;
  input.wheel_omega = nan;
  
  OUTPUT = zeros(rows,25);
  
  for ii=1:rows

    input.Fz    = INPUT(ii,1);
    input.kappa = INPUT(ii,2);
    input.alpha = INPUT(ii,3);
    input.gamma = INPUT(ii,4);
    input.path_curvature = INPUT(ii,5);
    input.Vx    = INPUT(ii,6);
    
    if (cols >= 7)
      input.pressure = INPUT(ii,7);
      if (cols >= 8)
        input.wheel_omega = INPUT(ii,8);
      end
    end

    output = tire.compute_steady(input, options);

    OUTPUT(ii,1:3) = output.force;
    OUTPUT(ii,4:6) = output.moment;
    OUTPUT(ii,7)  = output.kappa;
    OUTPUT(ii,8)  = output.alpha;
    OUTPUT(ii,9)  = output.gamma;
    OUTPUT(ii,10) = output.path_curvature;
    OUTPUT(ii,11) = vec3_dot(output.W(:,1), output.tire_contact_vel);
    OUTPUT(ii,12) = input.pressure; % not an output
    OUTPUT(ii,13) = output.rolling_radius;
    OUTPUT(ii,14) = output.rho;
    OUTPUT(ii,15) = output.contact_size(1);
    OUTPUT(ii,16) = output.trail(1);
    OUTPUT(ii,17:18) = output.grip;
    OUTPUT(ii,19) = output.wheel_omega;
    OUTPUT(ii,20) = output.loaded_radius;
    OUTPUT(ii,21:22) = output.relaxation;
    OUTPUT(ii,23) = output.contact_size(2);
    OUTPUT(ii,24) = output.vertical_stiffness;
    OUTPUT(ii,25) = output.phi;
  end
  
%  compare_params = compare_tire.get_param_vector();
%  result = tire.save([], [], 1, compare_params);
end

function d = vec3_dot(a, b)
  d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
end

