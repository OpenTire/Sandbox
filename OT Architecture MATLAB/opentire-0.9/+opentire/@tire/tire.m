%% opentire tire model abstract base class
classdef tire < matlab.mixin.Copyable

properties(Access = public)
  
  is_valid = false;   % true if all parameters are valid
  right_tire = false; % true if right tire, false otherwise
  filename = '';      % filename
  model_version = 0;  % model version number
  mdi_header = ''     % Extra ADAMS header information extracted from .TIR file
  
end % properties

methods(Abstract)
  new_options(tire)
  
  new_steady_input(tire)
  compute_steady(tire)
  
  new_dynamic_input(tire)
  compute_dynamic(tire)
end

methods(Abstract, Access=protected)
  on_params_changed(tire)
end

methods(Access = public)
  % Declare methods defined in separate files
  OUTPUT = mfeval(tire, INPUT, options)
  [ok, header, params] = load(tire, filename)
  ok = save(tire, filename, header, format, compare_params)
  result = diff(tire, compare_tire)
  result = utility(tire, command)
  
%% disp() - display a tire
%
function disp(tire)
  if (isempty(tire.filename))
    disp([class(tire) ' with reference params']);
  else
    disp([class(tire) ' from ' tire.filename]);
  end
end

%% Get number of parameters in parameter vector
function count = get_param_count(t)
  count = t.PARAM_INFO.param_count;
end

%% Get parameters as vector
function values = get_param_vector(t)
  count = t.PARAM_INFO.param_count;
  values = t.get_params_by_index((1:count)');
end

%% Set parameters from param vector
function set_param_vector(t, params)
  t.set_params_by_index((1:length(params))', params);
  t.on_params_changed();
end

function values = get_params_by_index(tire, indices)
  assert(false); % not implemented
end

function set_params_by_index(tire, indices, values)
  assert(false); % not implemented
end

end % methods(public)

methods(Static)
%  init_param_info();
end

properties(Constant)

%PARAM_INFO = init_param_info();

end

end % end classdef
