% plot - produce plot of tire inputs vs. outputs in the current figure.
%
% Inputs:
% tire - tire object
% inputs: cell array with names of x values and the x values to use
% outputs: cell array with names of output values to plot
% nominal_input - nominal input for all inputs not specified by xvalues.
%                 Default is for a left tire with nominal Fz, speed, and pressure.
% options - tire options struct (tire.new_options() if empty).
% compare_index - If non-zero, then used for plotting comparisons:
% compare_index 0: Normal plotting: one chart per window, each graph gets a unique color, with first always blue.
% compare_index > 0: Comparison plotting.  For each comparison chart, you call mftire_plot with a non-zero
%                   compare_index value for each comparison.  For a given compare_index all charts will have
%                   Caller must create figure or subplot to use.  Each chart will be drawn in the current figure
%                   with the same color -- you must create a figure or subplot before calling.
%
% IMPORTANT: You are responsible for creating a figure or subplot before calling this function.
%
% You can have 1 or 2 output values.  The first is shown on the left axis. The second, if present, is on the right.
% You can have 2 or 3 sets of input values:
%
% 2 inputs:
% For each first input:
%    for each second input:
%        plot output value(s)
%
% 3 inputs:
% For each first input:
%    For each second input:
%       For each third input
%           plot output values(s)
%
% NOTE: 
% If options.use_adapted_SAE is true, then the signs of path_curvature are reversed when plotting.
%
% See the code for a list of available input and output channel names.

% part of mftire 1.1.0
function plot(tire, inputs, outputs, nominal_input, options, compare_index)
  if nargin < 6
    compare_index = 0;
  end
  if nargin < 5 || isempty(options)
    options = tire.new_options();
  end
  if nargin < 4 || isempty(nominal_input)
    nominal_input = tire.new_nominal_input(false);
  end
  
  % The color blue is reserved for the first graph produced.
  % Afterward the colors cycle
  colors = 'brgkmrgkmrgkmrgkmrgkm';
  if compare_index >= 1
    % all comparison plots are plotted with same color
    colors(:) = colors(compare_index);
    
    % if drawing more than one curve, color first one differently
    if length(inputs{1,2}) > 1
      first_colors = 'cmkcmkcmkcmkcmkcmkcmk';
      colors(1) = first_colors(compare_index);
    end
  end
  
  Fz = nominal_input.Fz;
  alpha = nominal_input.alpha;
  kappa = nominal_input.kappa;
  gamma = nominal_input.gamma;
  curvature = nominal_input.path_curvature;
  press = nominal_input.pressure;
  Vcx = nominal_input.Vx;

  hold on; grid on;
  
  if size(inputs,1) == 3
    n = length(inputs{3,2});
  else
    n = length(inputs{2,2});
  end
  
  INPUTS = repmat([Fz alpha kappa gamma curvature Vcx press],n,1);
  
  if size(inputs,1) == 3
    params1 = inputs{1,2};
    params2 = inputs{2,2};
    params3 = inputs{3,2};
    
    for i1=1:length(params1)
      xx = zeros(0,1);
      yy = zeros(0,1);
      for i2 = 1:length(params2)
        
        INPUTS(:,input_chan_info(inputs{1,1})) = params1(i1);
        INPUTS(:,input_chan_info(inputs{2,1})) = params2(i2);
        INPUTS(:,input_chan_info(inputs{3,1})) = params3(:); % linspace(params3(1),params3(end),n)';
        
        OUTPUT = tire.mfeval(INPUTS, options);
        
        color = colors(i1);
        [x,label] = get_output_data(OUTPUT, outputs{1}, options);
        h = xlabel(label);
        h.Interpreter = 'none'; % so '_' doesn't get converted to funny characters in chart labels
        [y,label] = get_output_data(OUTPUT, outputs{2}, options);
        h = ylabel(label);
        h.Interpreter = 'none';
        
        if length(outputs) > 2
          yyaxis left;
        end
        plot(x,y,['-' color]);
        xx = [xx;x];
        yy = [yy;y];
        
        if length(outputs) > 2
          yyaxis right;
          [x,label] = get_output_data(OUTPUT, outputs{3}, options);
          h = xlabel(label);
          h.Interpreter = 'none';
          [y,label] = get_output_data(OUTPUT, outputs{4}, options);
          h = ylabel(label);
          h.Interpreter = 'none';
          plot(x,y,['-.' color]);
        end
      end
      % plot bounds of ellipse
      k = convhull(xx,yy);
      plot(xx(k),yy(k),['-' color]);

      % show results as we go
      drawnow('update');
    
    end
  else % 2 params
    params1 = inputs{1,2};
    params2 = inputs{2,2};
    
    for i1=1:length(params1)
      INPUTS(:,input_chan_info(inputs{1,1})) = params1(i1);
      INPUTS(:,input_chan_info(inputs{2,1})) = params2(:);
      
      OUTPUT = tire.mfeval(INPUTS, options);
      
      color = colors(i1);
      [x,label] = get_output_data(OUTPUT, outputs{1}, options);
      h = xlabel(label);
      h.Interpreter = 'none';
      [y,label] = get_output_data(OUTPUT, outputs{2}, options);
      h = ylabel(label);
      h.Interpreter = 'none';
      if length(outputs) > 2
        yyaxis left;
      end
      plot(x,y,['-' color]);
      
      if length(outputs) > 2
        yyaxis right;
        [x,label] = get_output_data(OUTPUT, outputs{3}, options);
        h = xlabel(label);
        h.Interpreter = 'none';
        [y,label] = get_output_data(OUTPUT, outputs{4}, options);
        h = ylabel(label);
        h.Interpreter = 'none';
        plot(x,y,['-.' color]);
      end
      
      % show results as we go
      drawnow('update');
      
    end
  end
  
end

% get named output channel info
function [index, name, units, scale] = output_chan_info(name)
  
  persistent info;
  if isempty(info)
    % first value is index into mfeval() OUTPUT vector
    info = {
      1, 'Fx', 'kgf', 1/9.81
      2, 'Fy', 'kgf', 1/9.81
      3, 'Fz', 'kgf', 1/9.81
      4, 'Mx', 'Nm', 1
      5, 'My', 'Nm', 1
      6, 'Mz', 'Nm', 1
      7, 'kappa', '-', 1
      7, 'lon_slip', '-', 1
      8, 'alpha', 'deg', 180/pi
      8, 'lat_slip', 'deg', 180/pi
      9, 'gamma', 'deg', 180/pi
      9, 'camber', 'deg', 180/pi
      9, 'inclination', 'deg', 180/pi
      10, 'curvature', '1/m', 1
      10, 'path_curvature', '1/m', 1
      11, 'Vx', 'km/h', 60*60/1000
      12, 'pressure', 'psi', 1/6895
      13, 'rolling_radius', 'mm', 1000
      13, 'Re', 'mm', 1000
      14, 'rho', 'm', 1
      15, 'contact_length', 'mm', 1000
      16, 'trail', 'mm', 1000
      17, 'mu_x', '-', 1
      18, 'mu_y', '-', 1
      19, 'omega', 'rpm', 60/(2*pi)
      20, 'loaded_radius', 'mm', 1000
      20, 'Rl', 'mm', 1000
      21, 'sigma_x', 'm', 1
      22, 'sigma_y', 'm', 1
      23, 'contact_width', 'mm', 1000
      24, 'stiffness', 'kfg/mm', 1/9.81/1000
      25, 'phi', '1/m', 1

      % Special normalized channels
      % index mod 100 => data channel to use
      % floor(index / 100) => data channel to divide by
      1 + 100*3, 'Fx/Fz', '-', 1
      2 + 100*3, 'Fy/Fz', '-', 1
      
      4 + 100*2, 'Mx/Fy', 'm', 1
      4 + 100*3, 'Mx/Fz', 'm', 1
      
      5 + 100*1, 'My/Fx', 'm', 1
      5 + 100*3, 'My/Fz', 'm', 1
      
      6 + 100*1, 'Mz/Fx', 'm', 1
      6 + 100*2, 'Mz/Fy', 'm', 1
      };
  end
  
  index = 0;
  units = '';
  scale = 1;
  for i=1:size(info,1)
    if strcmpi(name,info{i,2})
      index = info{i,1};
      name  = info{i,2};
      units = info{i,3};
      scale = info{i,4};
      break;
    end
  end
end

% Get named input channel info (just a column index for now)
function index = input_chan_info(name)
  switch lower(name)
  case 'fz'
    index = 1;
  case { 'kappa', 'lon_slip' }
    index = 2;
  case { 'alpha', 'lat_slip' }
    index = 3;
  case { 'gamma', 'camber' }
    index = 4;
  case { 'curvature', 'path_curvature' }
    index = 5;
  case 'vx'
    index = 6;
  case 'pressure'
    index = 7;
  case 'omega'
    index = 8;
  otherwise
    fprintf('can''t find channel %s\n', name);
    index = 0;
  end
end

% Get named output channel data
function [x, label, name, units] = get_output_data(output, name, options)
  [index,name,units,scale] = output_chan_info(name);
  if index == 0
    x = [];
    name = sprintf('%s undefined', name);
    units = '';
    label = name;
    return;
  end
  
  out_index = mod(index, 100);
  x = output(:,out_index) * scale;
  
  % HACK: Negative curvatures are used when making adapted_SAE plots.
  % Flip the sign here 
  if (options.use_adapted_SAE && out_index == 10) % curvature
    x = -x;
  end
  label = name;
  if ~isempty(units)
    label = sprintf('%s (%s)', name, units);
  end
  % handle normalization
  if index >= 100
    out_index = floor(index / 100);
    if (options.use_adapted_SAE && out_index == 10) % curvature
      x = -x;
    end
    x = x ./ output(:,out_index);
  end
  x = x(:);
  
end

