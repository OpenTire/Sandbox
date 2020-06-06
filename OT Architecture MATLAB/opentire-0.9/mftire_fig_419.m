% mftire_fig_419.m - Reproduce Fig. 4.19 on page 189 of TVD3

% part of mftire 1.1.0
function mftire_fig_419(tire)
  close all;
  
  if nargin == 0
    tire = opentire.mftire();%'tires/reference_tire_tvd3.tir';
  end
  
  if ischar(tire)
    tire = opentire.mftire(tire);
  else
    % copy the tire, so modifications to the tire won't affect tire of caller
    tire = copy(tire);
  end
  
  % NOTE: data from Table 4.2, p. 190
  % These parameters are DIFFERENT than reference_tire_tvd3.tir of Appendix 3
  
  tire.FNOMIN = 3000;
  tire.UNLOADED_RADIUS = 0.300;
  tire.Q_RE0 = 1.0; % rolling radius = loaded radius
  
  % Increase camber range so we can use book parameters
  tire.CAMMIN = -12*pi/180;
  tire.CAMMAX = 12*pi/180;
  tire.ALPMIN = -25*pi/180;
  tire.ALPMAX = 25*pi/180;
  
  % These parameters are changed relative to reference tire
  tire.QDRP2 = -1.5;
  tire.PECP1 = 0;
  
  left_tire = false;
  nominal = tire.new_steady_input(left_tire);
  nominal.Vx = 10/3.6;  % these are parking lot maneuvers
  nominal.Fz = 3000;
  nominal.gamma = 0*pi/180;
  nominal.path_curvature = 0;
  
  options = tire.new_options();
  options.want_turn_slip = true;
  
  % To repro charts in the book
  % Note that input parameter is always 1/path_radius, not -1/path_radius as in the book.
  options.use_adapted_SAE = true;
  
  steps = 150;
  
  % Plot some graphs from Fig. 4.19 on page 189 of VTD3 book
  
  half_contact_length = -0.1; % A nominal half contact patch length, for normalized curvature
  
  figure('Name','Fig. 4.19, page 189 (adapted SAE)');
  
  % Fx vs kappa
  subplot(3,2,1);
  title('turn slip Fy vs Fx, gamma = 0.0 and 8.0');
  
  nominal.gamma = 0*pi/180;
  nominal.path_curvature = 0;
  
  x_values = {
    'alpha' [-4 0 4 8 20]*pi/180
    'kappa' linspace(-.6,.6,100)
    };
  y_values = {'Fx' 'Fy' };
  tire.plot(x_values, y_values, nominal, options, 1);
  
  nominal.gamma = 8*pi/180;
  nominal.path_curvature = 0;
  tire.plot(x_values, y_values, nominal, options, 2);
  
  subplot(3,2,2);
  title('turn slip Mz vs Fx, gamma = 0.0 and 8.0');
  
  y_values = {'Fx' 'Mz' };
  nominal.gamma = 0*pi/180;
  nominal.path_curvature = 0;
  tire.plot(x_values, y_values, nominal, options, 1);
  
  nominal.gamma = 8*pi/180;
  nominal.path_curvature = 0;
  tire.plot(x_values, y_values, nominal, options, 2);
  
  % Fx vs kappa
  % Compare to: Figure 4.19: Fy vs. -a/R and Mz vs. -a/R
  subplot(3,2,3);
  title('Fig. 4.19: turn slip Fy vs kappa @ alpha, gamma = 0');
  
  nominal.gamma = 0*pi/180;
  
  x_values = {
    'alpha' [-4 0 4 8]*pi/180
    'curvature' linspace(0.0, 0.7, steps)/half_contact_length
    };
  y_values = {'curvature' 'Fy' };
  tire.plot(x_values, y_values, nominal, options);
  
  subplot(3,2,4);
  title('Fig. 4.19: turn slip Mz vs kappa @ alpha, gamma = 0');
  y_values = {'curvature' 'Mz' };
  tire.plot(x_values, y_values, nominal, options);
  
  % Fy/Mz vs alpha
  % Compare to: Figure 4.19: Fy vs. alpha and Mz vs alpha
  subplot(3,2,5);
  title('Fig. 4.19: turn slip Fy vs alpha @ curvature, gamma = 0');
  
  nominal.gamma = 0*pi/180;
  
  % Note that mftire_make_plot() flips the sign of curvature if use_adapted_SAE is true,
  % so plots look right.
  
  x_values = {
    'curvature' [-0.15 -0.10 -0.05 0.00 0.05 0.10 0.15 0.20 0.30 0.40 0.50 0.60 0.70]./half_contact_length
    'alpha' linspace(-15,15,steps)*pi/180
    };
  y_values = {'alpha' 'Fy' };
  tire.plot(x_values, y_values, nominal, options);
  
  y_values = {'alpha' 'Mz' };
  subplot(3,2,6);
  title('Fig. 4.19: turn slip Mz vs alpha @ curvature, gamma = 0');
  tire.plot(x_values, y_values, nominal, options);
end
