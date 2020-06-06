function brush_tire_test()
  
  tire = opentire.brush_tire(0);
  
  close all;

  plot_deflections(tire);
  plot_tire(tire);
end

function plot_tire(tire)
  
  options = tire.new_options();
  options.log_deflections = false;
  
  input = tire.new_steady_input(false);
  
  function [Fx,Fy,Mz,out] = brush_model(tire, Fz, Vc, kappa, alpha, gamma, curvature)
    
    input.Fz = Fz;
    input.Vx = Vc;
    input.kappa = kappa;
    input.alpha = alpha;
    input.gamma = gamma;
    input.path_curvature = curvature;

    output = tire.compute_steady(input, options);

    Fx = output.force(1);
    Fy = output.force(2);
    Mz = output.moment(3);
    out.p = output.p;
    out.e = output.e;
    out.b = output.b;
  end
  
  % nominal conditions
  Fz = tire.Fz;
  Vc = 30;
  curvature = 0;
  gamma = 0;
  
  % Fy,Mz vs. alpha
  kappa = linspace(0,-.2,5);
  cv = length(kappa);
  alpha = linspace(-25,25,100)*pi/180;
  ch = length(alpha);
  
  Fx = zeros(cv,ch);
  Fy = zeros(cv,ch);
  Mz = zeros(cv,ch);
  for iv = 1:cv
    for ih = 1:ch
      [Fx(iv,ih),Fy(iv,ih),Mz(iv,ih)] = brush_model(tire, Fz, Vc, kappa(iv), alpha(ih), gamma, curvature);
    end
  end
  
  figure;
  subplot(2,2,1); hold on; grid on;
  for iv = 1:cv
    plot(alpha,Fy(iv,:),'-');
  end
  
  subplot(2,2,2); hold on; grid on;
  for iv = 1:cv
    plot(alpha,Mz(iv,:),'-');
  end
  
  % Fy,Mz vs. Fx
  alpha = [-2 0 2 4 6 8]*pi/180;
  cv = length(alpha);
  kappa = linspace(-1,1,100);
  ch = length(kappa);
  
  Fx = zeros(cv,ch);
  Fy = zeros(cv,ch);
  Mz = zeros(cv,ch);
  for iv = 1:cv
    for ih = 1:ch
      [Fx(iv,ih),Fy(iv,ih),Mz(iv,ih)] = brush_model(tire, Fz, Vc, kappa(ih), alpha(iv), gamma, curvature);
    end
  end
  
  subplot(2,2,3); hold on; grid on;
  for iv = 1:cv
    plot(Fx(iv,:),Fy(iv,:),'-');
  end
  
  subplot(2,2,4); hold on; grid on;
  for iv = 1:cv
    plot(Fx(iv,:),Mz(iv,:),'-');
  end
  
  % Fy,Mz vs. Fx
  alpha = [-2 0 2 4 6 8]*pi/180;
  cv = length(alpha);
  kappa = linspace(-1,1,100);
  ch = length(kappa);
  
  Fx = zeros(cv,ch);
  Fy = zeros(cv,ch);
  Mz = zeros(cv,ch);
  for iv = 1:cv
    for ih = 1:ch
      [Fx(iv,ih),Fy(iv,ih),Mz(iv,ih)] = brush_model(tire, Fz, Vc, kappa(ih), alpha(iv), gamma, curvature);
    end
  end
  
  subplot(2,2,3); hold on; grid on;
  for iv = 1:cv
    plot(Fx(iv,:),Fy(iv,:),'-');
  end
  
  subplot(2,2,4); hold on; grid on;
  for iv = 1:cv
    plot(Fx(iv,:),Mz(iv,:),'-');
  end
end

function plot_deflections(tire)
  
  options = tire.new_options();
  options.log_deflections = true;
  
  input = tire.new_steady_input(false);
  
  function [Fx,Fy,Mz,out] = brush_model(tire, Fz, Vc, kappa, alpha, gamma, curvature)
    
    input.Fz = Fz;
    input.Vx = Vc;
    input.kappa = kappa;
    input.alpha = alpha;
    input.gamma = gamma;
    input.path_curvature = curvature;

    output = tire.compute_steady(input, options);

    Fx = output.force(1);
    Fy = output.force(2);
    Mz = output.moment(3);
    out.p = output.p;
    out.e = output.e;
    out.b = output.b;
  end
  
  Fz    = 3000;
  Vc    = 30;
  kappa = 0.00;       % long. slip (from -1 to ....)
  alpha = 4*pi/180;   % slip angle in degrees
  gamma = 0*pi/180;   % camber angle in degrees
  curvature = 1/1;    % lots of turnslip
  
  n = tire.elements;
  nrow = tire.rows;

  [~,~,~,out] = brush_model(tire, Fz, Vc, kappa, alpha, gamma, curvature);
  
  figure;
  
  cplot = 1 + nrow;
  subplot(cplot,1,1);
  hold on; grid on;
  
  if nrow == 1
    y_rows = 0;
  else
    y_rows = linspace(-0.5,0.5,nrow) * tire.contact_size(2) * tire.effective_width_scale;
  end
  
  for irow = 1:nrow
    
    b = squeeze(out.b(:,irow,:));
    e = squeeze(out.e(:,irow,:));
    
    % show row baseline and leading edge
    points = [b(1,end) b(2,1)
              b(1,end) y_rows(irow)
              b(1,1)   y_rows(irow)
              b(1,1)   b(2,1)
              b(1,end) b(2,1)];
    plot(points(:,1), -points(:,2),'r-');
    
    % plot carcass deflection in green
    points = [b(1,:)'  b(2,:)'
              b(1,1)   b(2,1)];
    plot(points(:,1), -points(:,2),'.-g');
    
    for i=1:size(b,2)
      points = [b(:,i) b(:,i)+e(:,i)];
      plot(points(1,:), -points(2,:),'b.-');
    end
  end
  
  for irow=1:nrow
    subplot(cplot,1,1+irow);
    hold on; grid on;
    x = squeeze(out.b(1,irow,:));
    p = squeeze(out.p(:,irow,:));
    norm_p = sqrt(sum(p(1:2,:).^2,1));
    y = [p(1,:); p(2,:); norm_p; p(3,:).*tire.mu0];
    plot(x,y, '.-');
  end
end
