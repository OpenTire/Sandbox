% mftire_plot - Plot tire data for one or more tires

% part of mftire 1.1.0
function mftire_plot(varargin)
    tires = varargin;
    % plot using reference tire
    if nargin == 0
      ta = opentire.mftire(0);
      %tires = { ta };
      tb = copy(ta);
      tb.LMUX = 0.90;
      tb.LMUY = 0.90;
      tires = {ta, tb};
    end
    tire_count = length(tires);
    
    % if tires passed as filenames, load them
    for i=1:tire_count
      if ischar(tires{i})
        tires{i} = opentire.mftire(tires{i});
      end
    end

    close all;

    compare_index = 1;
    if (tire_count == 1)
      compare_index = 0;
    else
      figure;
    end
    
    left_tire = false;
    
    steps = 50;
    
    Fz = linspace(1000,9000,3);
    gamma = linspace(-6,6,3)*pi/180;
    
    alpha = (-12:2:12)*pi/180;
    alpha_hi = linspace(-12,12,steps)*pi/180;
    
    kappa = linspace(-.2,.2,3);
    kappa_hi = linspace(-1,1,steps);
    
    for it = 1:tire_count
      if compare_index ~= 0
        compare_index = it;
      end
      tire = tires{it};
      
      options = tire.new_options();
      options.want_turn_slip = false;
      options.use_adapted_SAE = true;
      
      nominal = tire.new_steady_input(left_tire);
    
      % Fy,Mz vs alpha, varying kappa
      x_values = { 'kappa' kappa
                   'alpha' alpha_hi
                  };
      y_values = {'alpha' 'Fy' 'alpha' 'Mz'};
      open_figure(compare_index, 1, 'Fy,Mz vs. alpha @ kappa');
      tire.plot(x_values, y_values, nominal, options, compare_index);

      % Fy,Mz vs. alpha at varying Fz
      x_values = { 'Fz' Fz
                   'alpha' alpha_hi
                  };
      y_values = {'alpha' 'Fy' 'alpha', 'Mz'};
      open_figure(compare_index, 2, 'Fy,Mz vs. alpha @ Fz');
      tire.plot(x_values, y_values, nominal, options, compare_index);

      % Fx,Mz vs. alpha at varying gamma
      x_values = { 'gamma' gamma
                   'alpha' alpha_hi
                  };
      y_values = {'alpha' 'Fy' 'alpha' 'Mz'};
      open_figure(compare_index, 3, 'Fx,My vs. alpha @ gamma');
      tire.plot(x_values, y_values, nominal, options, compare_index);
      
      % Fy vs Fx grip ellipse, at varying Fz
      x_values = { 'Fz' Fz
                    'alpha' alpha
                    'kappa' kappa_hi
                    };
      y_values = {'Fx/Fz' 'Fy/Fz'};
      open_figure(compare_index, 4, 'Grip Ellipse @ Fz');
      tire.plot(x_values, y_values, nominal, options, compare_index);
    end
end

function open_figure(compare_index, chart_index, chart_title)
  if false && compare_index == 0
    % one chart per window
    if ~isempty(chart_title)
      figure('Name', chart_title);
    else
      figure;
    end
  else
    % comparison: 4 charts per figure
    subplot(2,2,chart_index);
    if ~isempty(chart_title)
      h = title(chart_title);
      % make sure '_' is printed properly
      h.Interpreter = 'none';
    end
  end
end