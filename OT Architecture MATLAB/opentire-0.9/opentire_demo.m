function opentire_demo()
  
%% OpenTire demo script

%% Load a tire, make a copy, change some parameters
tire = opentire.mftire('tires/reference_tire_tvd3.tir'); % load tire object
tire2 = copy(tire); % make a copy
tire2.LMUX = .8; % Modify some parameters
tire2.LMUY = .8;
%tire2.save('my_tire.tir'); % save .tir file

%% Show changed parameters
tire.diff(tire2);

%% Do a steady state computation
options = tire.new_options();

input = tire.new_steady_input(true);   % show a right tire
input.Fz = 5000;
input.Vx = 100/3.6;
input.kappa = .05;
input.alpha = 2*pi/180; % pos slip: neg Fy (to right)
input.gamma = 2.5*pi/180; % positive: leaning right
input.path_curvature = -1/10; % cornering to the right

fprintf('\n\nA steady-state force and moment computation\nSteady state conditions:\n');
fprintf('Fz = %.6g, Vx = %.2f km/h, kappa %.5g, alpha %.2f deg, gamma %.2f deg, path radius %.2f m\n', ...
  input.Fz, input.Vx*3.6, input.kappa, input.alpha*180/pi, input.gamma*180/pi, 1/input.path_curvature);

output = tire.compute_steady(input, options);

fprintf('\nForces and moments:\n');
fprintf('Fx = %.5g, Fy = %.5g, Fz = %.5g\n', output.force);
fprintf('Mx = %.5g, My = %.5g, Mz = %.5g\n', output.moment);
fprintf('\nSlip stiffnesses:\n');
fprintf('dFx/dkappa = %9.3g, dFy/dalpha = %9.3g\n', output.slip_stiffness(1:2));
fprintf('dFy/dgamma = %9.5g, dMz/dphi   = %9.5g\n', output.slip_stiffness(3:4));
fprintf('Loaded radius = %.2f mm\n', output.loaded_radius * 1000);
fprintf('Vertical stiffness = %.2f kgf/mm\n', output.vertical_stiffness*(1/9.81)/1000);
fprintf('Contact patch dimensions = %.2f x %.2f mm\n', output.contact_size*1000);
fprintf('Relaxation lengths = %.2f, %.2f mm\n', output.relaxation * 1000);

%% Plot some graphs for first tire
wait_for_key('Hit any key to see some graphs...');
mftire_plot(tire); % plot some graphs

%% Plot a comparison graph
wait_for_key('Hit any key to plot tire comparison...');
mftire_plot(tire, tire2);  % plot a comparison

%% Some tests
wait_for_key('Hit any key to run some tests...');
mftire_test(tire); % Some simple tests, sensitivity analysis

%% Compare computations with MFeval 4.0
wait_for_key('Hit any key to compare mftire output with MFeval 4.0 output...');
mftire_mfeval(tire); % Compare output of mftire with MFeval 4.0

%% Some other goodies
fprintf('Some other scripts...\n\n');

help mftire_fig_419
help mftire_loaded_radius_plot

% remove the tire file we created
%delete('my_tire.tir');

end

function wait_for_key(prompt)
  fprintf('\n');
  input(prompt);
  fprintf('\n');
end