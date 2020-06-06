% mftire_test - Runs some opentire.mftire() tests, shows sensitivity analysis

% part of mftire 1.1.0
function mftire_test(tire)
  
  if nargin == 0
    tire = 'tires/reference_tire_tvd3.tir';
    %tire = 'tires/MagicFormula62_Parameters.tir';
  end
  if ischar(tire)
    tire = opentire.mftire(tire);
  end
  if (~tire.is_valid())
    fprintf('Invalid tire parameters\n');
    return;
  end

  options = tire.new_options();
  options.want_turn_slip = false;
  options.transient_model = 0;
  
  road = opentire.plane_road_model();
  
  left_tire = tire;
  lin = left_tire.new_steady_input(false);
  
  right_tire = copy(tire); % make a copy
  rin = right_tire.new_steady_input(true);
  
  Fz = 3000;
  Vx = 100/3.6;
  kappa = 3*.01;
  alpha = 2*pi/180; % pos slip: neg fy (to right)
  gamma = 2.5*pi/180; % positive: leaning right
  curvature = -1/30; % cornering right
  
  [lin,rin] = fill_steady_input(lin,rin, Fz, alpha, gamma, kappa, Vx, curvature);

  ok = true;
  lout = left_tire.compute_steady(lin, options);
  print_output(lout);
  compute_stats(left_tire, lin, options);
  
  % make sure dynamic computation gives same results
  % as steady state for both left and right
  lin_d = make_dynamic_input(tire, lin, lout);
  
  states = [];
  algebraic_loops = [];
  
  [lout_d, states_dot, algebraic_loops] = left_tire.compute_dynamic(lin_d, road, states, algebraic_loops, options);
  if (~compare_output(lout, lout_d))
    ok = false;
    print_output(lout_d);
    fprintf('**** left dynamic != left static output\n\n');
  end

% loop to time dynamic evaluations
%   passes = 50000;
%   tic;
%   for pass=1:passes
%     lout_d = left_tire.compute_steady(lin_d, opt);
%   end
%   dt = toc*1e6/passes;
%   fprintf('dynamic: %.3f usec/call\n', dt);
%   return;

% test close_algebraic_loop
%   algloop_in = lin_d;
%   algloop_in.force = [0;0;0];
%   algout = left_tire.compute_steady(algloop_in, options);
%   [algloop_in, algout] = left_tire.close_algebraic_loop(algloop_in, algout, options, 10);

  % convert from dynamic back to static
  lin_2 = make_steady_input(left_tire, lin_d, lout_d);
  lout_2 = left_tire.compute_steady(lin_2, options);
  
  if (~compare_output(lout, lout_2))
    ok = false;
    print_output(lout_2);
    fprintf('**** left round trip static != left static output\n\n');
  end
  
  rout = right_tire.compute_steady(rin, options);
  print_output(rout);
  compute_stats(right_tire, rin, options);

  rin_d = make_dynamic_input(right_tire, rin, rout);
  
  [rout_d, states_dot, algebraic_loops] = right_tire.compute_dynamic(rin_d, road, states, algebraic_loops, options);
  
  if (~compare_output(rout, rout_d))
    ok = false;
    print_output(rout_d);
    fprintf('**** right dynamic != right static output\n\n');
  end
  
  rin_2 = make_steady_input(right_tire, rin_d, rout_d);
  rout_2 = right_tire.compute_steady(rin_2, options);
  if (~compare_output(rout, rout_2))
    ok = false;
    print_output(rout_2);
    fprintf('**** right round trip static != right static output\n\n');
  end
  
  lout_r = swap_left_right_output(rout_2);
  
  % Make sure right output == left output when reflected about YZ plane
  if (~compare_output(lout, lout_r))
    ok = false;
    print_output(lout_r);
    fprintf('**** right output != reflected left output\n\n');
  end
  
  if ~ok
    fprintf('\nTest FAILED\n');
  else
    fprintf('\nTest successful\n');
  end
  
end

function swapped = swap_left_right_output(output)
  swapped = output;
  
  swapped.right_tire = ~swapped.right_tire;
  
  swapped.force(2)  = -swapped.force(2);
  swapped.moment(1) = -swapped.moment(1);
  swapped.moment(3) = -swapped.moment(3);
  
  swapped.alpha            = -swapped.alpha;
  %swapped.transient_dot(1) = -swapped.transient_dot(1);
  swapped.tire_slip_vel(2) = -swapped.tire_slip_vel(2);

  swapped.gamma = -swapped.gamma;
  swapped.path_curvature = -swapped.path_curvature;
  
  swapped.R = diag([1 -1 1]) * swapped.R * diag([1 -1 1]);
  swapped.W = diag([1 -1 1]) * swapped.W * diag([1 -1 1]);

  swapped.tire_contact(2) = -swapped.tire_contact(2);
  swapped.trail(2) = -swapped.trail(2);
end


function ok = compute_stats(tire, in, opt)
  ok = false;
  psi_to_pa = 6895;
  pa_to_psi = 1/psi_to_pa;
  
  x0 = [in.Fz; in.Vx; in.kappa; in.alpha; in.gamma; in.path_curvature; in.pressure];
  
  function [y,out] = do_eval(x)
    in.Fz = x(1);
    in.Vx = x(2);
    in.kappa = x(3);
    in.alpha = x(4);
    in.gamma = x(5);
    in.path_curvature = x(6);
    in.pressure = x(7);
    out = tire.compute_steady(in, opt);
    y = [out.force; out.moment; out.loaded_radius];
  end
  
  [y0, out] = do_eval(x0);
  
  fprintf('\n\nSteady State\n');
  fprintf('Fz = %.0f N, Vx = %.1f km/h, kappa = %.1f%%, alpha = %.1f deg, gamma = %.1f deg, curvature = %.3f, press = %.1f psi\n\n', ...
    x0 .* [1; 3.6; 100; 180/pi; 180/pi; 1; 1/psi_to_pa]);
  
  J = two_sided_jacobian(@do_eval, x0, y0);
  
  % convert units of denominator
  J = J * diag(1 ./ [1/1000 3.6/10 100 180/pi 180/pi 1/.1 pa_to_psi]);
  % convert loaded_radius units
  J(7,:) = J(7,:) * 1000; % m -> mm
  
  fprintf('\n');
  fprintf('       SENSITIVITIES:       Fx       Fy       Fz       Mx       My       Mz   ld_radius\n');
  fprintf('    Fz effect per kN: %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n', J(:,1));
  fprintf('Vx effect per 10 kmh: %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n', J(:,2));
  fprintf('  Kappa effect per %%: %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n',J(:,3));
  fprintf('Alpha effect per deg: %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n', J(:,4));
  fprintf('Gamma effect per deg: %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n', J(:,5));
  fprintf('Phi effect per 1/.1m: %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n', J(:,6));
  fprintf('Pressure eff per psi: %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n', J(:,7));
end

function [lin,rin] = fill_steady_input(lin,rin, Fz, alpha, gamma, kappa, Vx, curvature)
  lin.Fz = Fz;       rin.Fz = Fz;
  lin.Vx = Vx;       rin.Vx = Vx;
  lin.kappa = kappa; rin.kappa = kappa;
  lin.alpha = alpha; rin.alpha = -alpha;
  lin.gamma = gamma; rin.gamma = -gamma;
  lin.loaded_radius = nan; rin.loaded_radius = nan;
  lin.path_curvature = curvature; rin.path_curvature = -curvature;
end

% Create a new input, from from_input and its output from_output,
% converting from steady-state to dynamic or vice versa.
% The new input should give the same output as the original from_input.
function input = make_dynamic_input(tire, from_input, from_output)
  input = tire.new_dynamic_input;
  
  input.pressure = from_input.pressure;
  input.right_tire = from_input.right_tire;
  input.grip_adjust = [1;1];
  
  % input was steady state:
  % convert to dynamic input
  
  Vx = dot(from_output.tire_contact_vel, from_output.W(:,1));
  input.frame_omega = [0;0;Vx * from_output.path_curvature];

  input.wheel_omega = from_output.wheel_omega;

  % Compute wheel center position and velocity
  radial = from_output.R(:,3)*-from_output.loaded_radius; % vector from center to contact
  input.wheel_center = from_output.tire_contact - radial;

  % correct wheel center velocity for frame rotation
  input.wheel_center_vel = from_output.tire_contact_vel - cross(input.frame_omega, radial);

  input.wheel_axis = from_output.R(:,2);
end

% Create a new steady input, from from_input and its output from_output,
% The new input should give the same output as the original from_input.
function input = make_steady_input(tire, from_input, from_output)
  input = tire.new_steady_input;
  
  input.right_tire = from_input.right_tire;
  input.pressure = from_input.pressure;
  input.Fz = from_output.force(3);
  input.Vx = dot(from_output.tire_contact_vel, from_output.W(:,1));
  input.kappa = from_output.kappa;
  input.alpha = from_output.alpha;
  input.gamma = from_output.gamma;
  input.path_curvature = from_output.path_curvature;
  input.loaded_radius = nan;
end

function print_output(out)
  psi_to_pa = 6895;
  pa_to_psi = 1/psi_to_pa;
  
  side = 'LEFT'; if out.right_tire, side = 'RIGHT'; end
  fprintf('\n');
  fprintf('%20s: %s\n', 'TIRE', side);
  wheel_center = out.tire_contact + out.R(:,3)*out.loaded_radius;
  fprintf('\n');
  fprintf('%20s: %8.2f %8.2f %8.2f N\n', 'force', out.force);
  fprintf('%20s: %8.2f %8.2f %8.2f N/m\n', 'moment', out.moment);
  fprintf('\n');
  
  fprintf('%20s: %6.2f %%\n', 'kappa', out.kappa*100);
  fprintf('%20s: %6.2f deg\n', 'alpha', out.alpha*180/pi);
  fprintf('%20s: %6.2f deg\n', 'gamma', out.gamma*180/pi);
  fprintf('%20s: %7.3f 1/m\n', 'phi', out.phi);
  fprintf('%20s: %7.3f 1/m\n', 'path_curvature', out.path_curvature);
  path_radius = 0; if isfinite(out.path_curvature) && out.path_curvature ~= 0, path_radius = 1/out.path_curvature; end
  fprintf('%20s: %6.2f m\n', 'path_radius', path_radius);
  %fprintf('%20s: %5.1f psi\n', 'pressure', out.pressure * pa_to_psi);
  fprintf('%20s: %6.2f rad/s\n', 'wheel_omega', out.wheel_omega);
  fprintf('%20s: %6.2f mm\n', 'rolling_radius', out.rolling_radius*1000);
  fprintf('%20s: %6.2f mm\n', 'loaded_radius', out.loaded_radius*1000);
  %fprintf('%20s: %6.2f psi\n', 'pressure', out.pressure*pa_to_psi);
  fprintf('\n');
  fprintf('%20s: %6.2f mm\n', 'rho', out.rho*1000);
  fprintf('%20s: %6.2f kgf/m\n', 'stiffness', out.vertical_stiffness/9.81/1000);
  fprintf('%20s: %6.2f %6.2f mm\n', 'trail', out.trail*1000);
  fprintf('%20s: %6.2f %6.2f\n', 'grip', out.grip);
  fprintf('%20s: %5.1f %5.1f  mm\n', 'relaxation', out.relaxation*1000); %./(2*pi*out.rolling_radius)*100);
  fprintf('%20s: %5.1f %5.1f %6.1f %6.1f (N/%%, N/deg, N/deg, .1)\n', 'slip_stiffness', out.slip_stiffness .* [1/100; 1/(180/pi); 1/(180/pi); .1]);
  fprintf('%20s: %6.2f %6.2f mm\n', 'contact_size', out.contact_size*1000);
  fprintf('%20s: %5.1f %5.1f %6.1f\n', 'wheel center', wheel_center*1000);
  fprintf('%20s: %5.1f %5.1f %6.1f\n', 'tire contact', out.tire_contact*1000);
  fprintf('%20s: %5.1f %5.1f %6.1f m/s\n', 'tire contact vel', out.tire_contact_vel);
  fprintf('%20s: %5.1f %5.1f m/s\n', 'tire slip vel', out.tire_slip_vel);
  %fprintf('%20s: %6.4g %6.4g\n', 'transient_dot', out.transient_dot);
end

function ok = compare(name, out1, out2)
  v1 = out1.(name);
  v2 = out2.(name);
  
  rel_tol = 1e-5;
  abs_tol = 1e-4;
  
  ok = true;
  if isscalar(v1)
    d = abs(v1-v2) / max(abs_tol, max(abs(v1),abs(v2)));
    if d >= rel_tol
      fprintf('%s: VALUES DIFFERENT:\n', name);
      fprintf('    1: %9.6g\n', v1);
      fprintf('    2: %9.6g\n', v2);
      fprintf(' diff: %9.3g (%.3g)\n', v2-v1, d);
      ok = false;
    end
  else
    d = norm(v1-v2) / max(abs_tol, norm(max(abs(v1),abs(v2))));
    if d >= rel_tol
      fprintf('%s: VALUES DIFFERENT:\n    1: ', name);
      fprintf('%9.6g ', v1(:)');
      fprintf('\n    2: ');
      fprintf('%9.6g ', v2(:)');
      fprintf('\n diff: ');
      fprintf('%9.3g ', v2(:)'-v1(:)');
      fprintf(' (%.3g)\n', d);
      ok = false;
    end
  end
end

function ok = compare_output(out1, out2)
  ok = true;
  f = fields(out1);
  
  for i=1:length(f)
    ok = ok & compare(f{i}, out1, out2);
  end
end

function J = two_sided_jacobian(f, x, y)
  m = length(y);
  n = length(x);
  J = zeros(m, n);
  delta = 64*sqrt(eps);
  for j=1:n
    xj = x(j);
    h = max(delta, delta * abs(xj));
    x(j) = xj + h;
    xp = x(j);
    yp = f(x);
    
    x(j) = xj - h;
    ym = f(x);
    
    J(:,j) = (yp - ym) ./ (xp - x(j));
    
    x(j) = xj;
  end
end
