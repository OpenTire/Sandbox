%% diff - Compare two tire objects
% Prints differing parameters to console window.
%
% part of mftire 1.1.0
function result = diff(tire, compare_tire)
  compare_params = compare_tire.get_param_vector();
  result = tire.save([], [], 0, compare_params);
end
