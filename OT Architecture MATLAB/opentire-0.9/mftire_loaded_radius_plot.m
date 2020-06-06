% mftire_loaded_radius_plot - Shows MF 6.2 loaded radius effects of alpha, kappa, camber
% Plots Fz at various conditions with fixed loaded radius.

% part of mftire 1.1.0
function mftire_loaded_radius_plot(tire)
% Plots Fz at various conditions with fixed loaded radius.
  
  if nargin == 0
    tire = 'tires/MagicFormula62_Parameters.tir';
  end
  if ischar(tire)
    tire = opentire.mftire(tire);
  end

  options = tire.new_options();
  
  input = tire.new_steady_input();
  
  Fz = [0 1000 3000];% 5000 7000 9000];
  alpha = linspace(-5,5,50)*pi/180;
  kappa = linspace(-.2,.2,50);
  gamma = linspace(-1,1,50)*pi/180;
  
  Fz_a = zeros(length(Fz),length(alpha));
  input.kappa = 0;
  for ifz = 1:length(Fz)
    input.Fz = Fz(ifz);
    input.alpha = 0;
    input.kappa = 0;
    input.gamma = 0;
    input.loaded_radius = nan;
    input.loaded_radius = tire.compute_steady(input, options).loaded_radius;
    
    for ial = 1:length(alpha)
      input.alpha = alpha(ial);
      output = tire.compute_steady(input, options);
      Fz_a(ifz,ial) = output.force(3);
    end
  end
  
  Fz_k = zeros(length(Fz),length(alpha));
  input.alpha = 0;
  for ifz = 1:length(Fz)
    input.Fz = Fz(ifz);
    input.alpha = 0;
    input.kappa = 0;
    input.gamma = 0;
    input.loaded_radius = nan;
    input.loaded_radius = tire.compute_steady(input, options).loaded_radius;
    
    for ik = 1:length(kappa)
      input.kappa = kappa(ik);
      output = tire.compute_steady(input, options);
      Fz_k(ifz,ik) = output.force(3);
    end
  end
  
  Fz_g = zeros(length(Fz),length(gamma));
  input.alpha = 0;
  for ifz = 1:length(Fz)
    input.Fz = Fz(ifz);
    input.alpha = 0;
    input.kappa = 0;
    input.gamma = 0;
    input.loaded_radius = nan;
    input.loaded_radius = tire.compute_steady(input, options).loaded_radius;
    
    for ig = 1:length(gamma)
      input.gamma = gamma(ig);
      output = tire.compute_steady(input, options);
      Fz_g(ifz,ig) = output.force(3);
    end
  end
  
  close all;
  figure('Name','Constant Loaded Radius');
  
  subplot(3,1,1); hold on; grid on;
  plot(alpha*180/pi,Fz_a,'.-');
  title('Fz vs alpha @ constant radius at 0,1,3,5,7,9 kN');
  
  subplot(3,1,2); hold on; grid on;
  plot(kappa*100,Fz_k,'.-');
  title('Fz vs kappa @ constant radius at 0,1,3,5,7,9 kN');
  
  subplot(3,1,3); hold on; grid on;
  plot(gamma*180/pi,Fz_g,'.-');
  title('Fz vs camber @ constant radius at 0,1,3,5,7,9 kN');
end
