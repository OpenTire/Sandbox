% tyr501Harty_2013.m
% (c) Damian Harty 24 oct 1999-2020
%
% Translated from tyr501Harty_2013.f
% NK - 29-May-2020 - translation version 2.0
%
% Works completely in ISO/Tydex-W. Road model and tire kinematics
% code taken from OpenTire prototype.
%
% Major changes over 2013 version:
% 1) Variable names were changed to more closely match my other tire model code (mftire)
%    in preparation for dropping this model into opentire.  My apologies if these
%    changes render this code unreadable to you.
% 2) The unloaded radius is now used for the effective rolling radius;
%    the 2013 model used the loaded radius.
% 3) I discovered that very high camber angles had big problems
%    in my first translation: clearly I broke something.
%    However I was unable to grok the old code and no amount
%    of sign flipping gave satisfactory results.  So I had a
%    play and came up with a model that applies Fy linearly with
%    gamma up to a breakaway point, then smoothly transitions back to
%    zero at a settable rate.  I tried to match the gains, etc, of the old code,
%    at least in the region abs(gamma) < peak gamma.
% 4) I also poked around in the friction ellipse code, trying to make
%    transitions a little smoother, etc.  In the end I believe this
%    code works like the original code.  I was trying to adjust things
%    so that the peak mu values were actually reached in both Fx and Fy,
%    to get rid of the clipping seen as flat edges of the grip ellipse.
%
% ProDrive Concept Tyre Model
%
% A quick & dirty tyre model originally designed to plug into ADAMS
% as the Fiala model does.
%
% Unlike Fiala, critical slip angle is broadly independent
% of load and initial cornering stiffness is strongly
% load dependent.
%
% These attributes better represent a modern radial tyre
% than does either the Fiala or University of Arizona model.
%
% The model does handle comprehensive slip.  lateral force
% generation is zero at peak longitudinal force slip ratio
% (typically about 20%) but returns to a value around one
% tenth of the peak lateral force as the wheel progresses
% beyond that limit.  this may result in poor post-spin
% performance.  the force generated with locked wheels is
% aligned with the wheel plane; this is incorrect.
%
% Longitudinal force generation is assumed to be symmetric
% for tractive and braking slip.  this is not generally
% true beyond the critical slip ratio for real tyres but
% is reasonable up to that point.  this tyre will over
% estimate longitudinal forces for tractive slip and
% slightly underestimate them for braking slip in the
% post-critical regions.
%
% Camber thrust is included as for the motorcycle tire
% model using "taut string" logic.  lateral migration of
%NK: Is the model in use still a "taut string" model?
% the contact patch is also included, as for the motorcycle
% tyre model.
%
% Aligning torque calculation includes the lateral force
% due to camber.  this is not quite right as the camber
% force mechanism has no pneumatic trail associated with
% it.  pay attention if using this for motorcycle work;
% consider reworking it so that Mz does not include the
% camber force.  the form of the aligning torque is a
% bit poor and would benefit from some more thought;
% pneumatic trail collapses linearly with lateral force.
%
function [wheel_center_force, wheel_center_moment, ...
    force, moment, ...
    states_dot, ...
    error_code] ...
    = tyr501Harty_2020(tire, ...
    computation_code, solver_status, ...
    time, ...
    frame_omega, ...
    wheel_center, T, wheel_center_vel, ...
    wheel_angle, wheel_omega, ...
    states)

if (nargin == 0), test(); return; end

% imode values
dynamic = 1;
static  = 0;

% jobflg values from the integrator/solver passed through the Standard Tire Interface (STI).
% jobflg = 0 = normal   Normal Mode
% jobflg = 1 = inquire  Subroutine must return the actually used dimensions of workareas
% jobflg = 2 = init     First initialization before starting the simulation
% jobflg = 3 = reset    Re-initialization during simulation
% jobflg = 4 = sstep    Successful step (not used in A/Car)
% jobflg = 5 = diff     DFLAG is true (doing numeric differencing) (A/Car deviation)
% jobflg = 99 = endsim  Final call of tire model (not used in A/Car)

noeval  = -1;
normal  = 0;
inquire = 1;
init    =  2;
reset   =  3;
sstep   =  4;
diff    =  5;
endsim  =  99;

% error code
error_code = 0;

% Get model parameters
use_mode      = tire.USE_MODE;

% General tire parameters
Fz0           = tire.RATED_LOAD;

R_0           = tire.UNLOADED_RADIUS;
tire_width    = tire.WIDTH;
vertical_stiffness    = tire.VERTICAL_STIFFNESS;
vertical_damping      = tire.VERTICAL_DAMPING;

% Mx model parameters
CMyFz         = tire.ROLLING_RESISTANCE;

C_alpha       = tire.CALPHA; % called CMX in .TIR file
C_gamma       = tire.CGAMMA;
mu_min        = tire.UMIN;
mu_max        = tire.UMAX;
relaxation_length = tire.RELAXATION_LENGTH;
peak_alpha    = tire.ALPHA_CRITICAL * pi/180; % convert from degrees
ay            = tire.CURVATURE_FACTOR_ANGLE;

mu_y_scale0   = tire.SCALE_FACTOR_LATERAL;
CMuyFz        = tire.SCALE_FACTOR_DIM;

peak_kappa    = tire.SLIP_RATIO_CRITICAL * 0.01; % convert from percent
ax            = tire.CURVATURE_FACTOR_RATIO;
trail_scale   = tire.PNEUM_TRAILING_SCALING;
camber_Fy_moment_arm = tire.PNEUMATIC_LEAD_CAMBER;
peak_gamma_threshold = tire.LIMIT_CAMBER_ONSET_FRIC;

%   if ( solver_status == init || solver_status == reset )
%   end

% initialize mode (static or dynamic)
imode = dynamic;
if ( computation_code == 0 ), imode = static; end

% set flag for quasi-static analyses
staflg = false;
if ( computation_code == 2 ), staflg = true; end

% startup smoothing function:
%
% the mdi tire models include a feature for smoothing the
% tire forces around time = 0.0.  so, for example, if there's
% some initial slip angle at time=0.0, the lateral force
% builds up slowly instead of acting like a step input.
% this helps the integrator get started.

%   startup_scale = 1.0;
%   if (use_mode >= 2 && ~staflg)
%     end_startup_scale_time = 0.5;
%     startup_scale = step(time, 0.0,0.0, end_startup_scale_time,1.0);
%   end

rolling_radius = R_0;

% compute tire contact point and other aspects of road contact
[rho, tire_contact, road_normal, road_grip] ...
    = road_model(solver_status, time, T, wheel_center, R_0, tire_width);

% If we're not doing normal simulation or finite-differencing,
% treat as if the wheel has come off the ground.
if (~(solver_status == normal || solver_status == diff))
    rho = 0;
end

[Vc, Vs, W, kappa, alpha, gamma] = tire_kinematics(frame_omega, T, wheel_center_vel, wheel_omega, rolling_radius, road_normal);

% Tire relaxation model
% lag the slips for tire relaxation effects:

% If the wheel has come off the ground, the effective slip angles are zero.
% It would be nice to reset the states to zero when the tire comes off
% the ground, so that the relaxation effect restarts from zero when
% the tire touches ground again, but alas these are managed by the
% calling solver and can only be affected via states_dot.
if (rho <= 0)
    kappa = 0;
    alpha = 0;
    states_dot = [0;0];
    
    Vs = [0;0];  % Reset slip velocities too.
    
else
    
    kappa_unlagged = kappa;
    alpha_unlagged = alpha;
    states_dot = [0;0]; % default to zero
    
    tiny_relaxation = 1e-4;
    if (relaxation_length > tiny_relaxation && imode ~= static )
        gain  = abs(Vc(1)) / relaxation_length;
        kappa = states(1);
        alpha = states(2);
        states_dot(1) = gain * (kappa_unlagged - kappa);
        states_dot(2) = gain * (alpha_unlagged - alpha);
    end
end

tan_alpha = tan(alpha);

% normalized slip quantities
nkappa = abs(kappa / peak_kappa);
nalpha = abs(alpha / peak_alpha);

% calculate normal loads due to stiffness (always <= zero) --
% loaded_radius = R_0 - rho;
Fz_spring = rho * vertical_stiffness;

% calculate normal loads due to damping --
Fz_damping = Vc(3) * vertical_damping;

Fz = max(Fz_spring + Fz_damping, 0.0);

if (imode == static || rho == 0)
    
    % Solving for vertical steady state, or wheel is off the ground
    % zero everything but fz (which will be zero if rho == 0)
    Fx = 0; Fy = 0;
    Mx = 0; My = 0; Mz = 0;
    
else % (imode == dynamic)
    
    % grip rolloff as a function of combined slip.
    % (chosen to give similar results to original
    % definition of sqrt(kappa^2 + tan_alpha^2)
    
    % allows changing proportion of alpha in
    % mu_per_combined_slip computation
    mu_tan_alpha_scale = 1;
    
    % combined mu_slip
    mu_per_combined_slip = sqrt( kappa^2 + (tan_alpha*mu_tan_alpha_scale)^2 );
    
    mu = mu_max + (mu_min - mu_max) * mu_per_combined_slip;
    
    % modify coefficient of friction based on road surface factor:
    mu = mu * road_grip;
    mu_x = mu;
    
    % Lateral grip scale, diminished with Fz
    mu_y_scale = mu_y_scale0 + (Fz - Fz0) * CMuyFz;
    
    mu_y = mu * mu_y_scale;
    
    % max Fx and Fy
    Fx_max = Fz * mu_x;
    Fy_max = Fz * mu_y;
    
    % Longitudinal force
    
    if (nkappa <= 1)
        % -- exponential rise (1-exp(-ax * nkappa)) below peak_kappa --
    else
        % -- clamped at sliding friction above.
        % This results in a linear decay due to the
        % combined-slip dependent reduction in the calculation of mu.
        nkappa = 1;
    end
    
    Fx = Fx_max * SGN(kappa) * (1 - exp(-ax * nkappa));
    
    % Lateral force
    Fy = 0;
    
    tiny_alpha = 1e-8;
    if (abs(alpha) > tiny_alpha)
        
        nalpha = abs(alpha / peak_alpha);
        
        if (abs(nalpha) <= 1)
            % -- exponential rise (1-exp(-ay)) below critical slip angle --
        else
            nalpha = 1;
        end
        % Positive alpha produces negative Fy
        Fy = Fy_max * -SGN(alpha) * (1 - exp(-ay * nalpha));
    end
    
    % aligning torque based on intermediate fy excluding camber force.
    
    % -- contact patch length --
    contact_length = sqrt(R_0^2 - (R_0 - rho)^2) * 2;
    
    % -- trail_patch_size_factor is because lever arm is not the entire contact patch length. --
    
    % -- parameter trail_scale should be set to 1.0 for tyres with rectangular
    %    contact patches (i.e. car tyres) and 0.5 for tyres with elliptical
    %    contact patches (i.e. motorcycle tyres.) --
    
    base_trail_factor = 1/6;
    trail_moment_arm = base_trail_factor * trail_scale;
    trail = contact_length * trail_moment_arm * max(0, 1 - nalpha);
    
    % positive trail is behind the tire contact center.
    % So positive trail and positive Fy produce a negative Mz.
    Mz = Fy * -trail;
    
    % -- camber_breakaway represents aggression of departure at limit; high value
    %    implies high limit & aggressive departure, lower value implies
    %    progression.
    
    % Camber thrust effects. Additional Fy is provided by camber
    % linearly until peak_mu_camber
    % (old code increased as function of tan(camber), nearly the same thing)
    
    peak_mu_camber = mu * peak_gamma_threshold /C_gamma;
    peak_gamma = atan(peak_mu_camber);
    
    %peak_gamma = 30*pi/180;
    camber_breakaway = 15*pi/180; % camber breakaway when peak_camber < camber > peak_camber + camber_breakaway
    
    abs_gamma = abs(gamma);
    if (abs_gamma < peak_gamma)
        s = abs_gamma / peak_gamma;
        s = smoothstep(0.5 + s/2)*2 - 1;
    else
        s = (abs_gamma - peak_gamma) / camber_breakaway;
        s = 1 - min(s, 1);
        s = smootherstep(s); % sharper but soft at both ends
    end
    mu_camber = s * peak_mu_camber;
    
    % Positive camber and positive Fz produce negative Fy.
    Fy_camber = mu_camber * Fz * -SGN(gamma);
    
    Fy = Fy + Fy_camber;
    
    Fy_pure = Fy; % save pure fy
    
    % mitigate fy depending on "friction ellipse"
    sign_Fx = SGN(Fx);
    sign_Fy = SGN(Fy);
    
    % normalized friction
    nFx = Fx/Fx_max;
    nFy = Fy/Fy_max;
    
    Fr = nFx^2 + nFy^2;
    if (Fr > 1.0)
        F_scale = 1 / sqrt(Fr);
        Fx1 = Fx * F_scale;
        Fy1 = Fy * F_scale;
        
        % -- Alternate formulation used with higher slip ratios that
        % reduces the maximum revised formulation for highest slip ratios arrived at
        % by consideration of contact patch velocity.  gives pleasing
        % results for wheels locked and wheels spinning cases.
        
        % used to change proportion of alpha vs. kappa speed effect
        speed_effect_alpha_scale = 1;
        speed_effect_tan_alpha = tan_alpha * speed_effect_alpha_scale;
        
        speed_effect = 1 / (1 + (speed_effect_tan_alpha / kappa)^2);
        %speed_effect = 1 / (1 + (tan_alpha / kappa)^2);
        
        % Basically this formula reduces Fx2 to a speed-dependent
        % maximum, and Fy2 gets whatever's left over in a Pythagorean sense.
        Fx2 = sqrt( Fx_max^2 * speed_effect ) * sign_Fx;
        Fy2 = sqrt( Fx_max^2 - Fx2^2 ) * sign_Fy;
        
        % Blend our two friction ellipse functions.
        
        % For kappa below kappa_lo, use [Fx1,Fy1],
        % for kappa above kappa_hi, use [Fx2,Fy2].
        % Between kappa_lo and kappa_hi, interpolate between
        % Fx1,Fy1 and Fx2,Fy2 according to the step() function.
        kappa_lo = 0.50;
        kappa_hi = 1.00;
        
        interpolate = step(abs(kappa), kappa_lo,0.0, kappa_hi,1.0);
        
        Fx = Fx1 + (Fx2 - Fx1) * interpolate;
        Fy = Fy1 + (Fy2 - Fy1) * interpolate;
        
        % -- mitigate camber forces too, for subsequent aligning moment calculations
        
        %     is this right? doesn't it significantly corrupt aligning torque for
        %     a locking wheel at a high slip angle?
        if (Fy_pure ~= 0) % avoid division by zero
            Fy_camber = Fy_camber * Fy/Fy_pure;
        end
        
    end % if (fr_ellip > 1.0)
    
    % -- reverse sign of Fx if we are moving backwards --
    Fx = Fx * SGN(Vc(1));
    
    % rolling resistance moment (same as Fiala tire)
    My = CMyFz * -Fz;
    
    % Compute righting moment due to lateral contact patch shift (mx)
    % Use ca as "shape factor" to add to or subtract righting moment
    % from ADAMS' toroidal assumption.  ca > 1 = fatter than toroid
    % ca < 1 = more like blade.
    
    % add aligning torque based on lateral offset of contact patch and
    % longitudinal forces to give "stand up under braking" behaviour for
    % motorcycles or tramlining for cars.
    
    % Positive gamma will produce a negative offset,
    % which will produce a negative Mx with a positive Fz,
    % and a positive Mx with positive driving Fx.
    
    htw = tire_width / 2; % half tread width;
    camber_FzFx_moment_arm = -2 * gamma * htw * (C_alpha - 1);
    
    % don't allow offset past physical edge of tire
    if (abs(camber_FzFx_moment_arm) > htw)
        camber_FzFx_moment_arm = htw * SGN(camber_FzFx_moment_arm);
    end
    
    Mx = Fz * camber_FzFx_moment_arm;
    
    Mz = Mz - Fx * camber_FzFx_moment_arm;
    
    % Measured data shows evidence of significant "pneumatic lead" on
    % camber force data, aligning moment further modified to reflect this.
    % real data shows small dependency on load, some dependency on camber
    % angle at low cambers; constant lead formulation neglects load
    % dependency and may overestimate torques at small cambers. however,
    % camber forces are low and so torques are low too.
    
    % A pneumatic "lead" is a positive offset of Fy_camber moment arm,
    % ahead of the tire contact point, which will produce positive Mz
    % for positive values of Fy.
    Mz = Mz + Fy_camber * camber_Fy_moment_arm;
    
end % end of imode == dynamic

% Forces and moments at contact point, in the Tydex W frame
force  = [Fx;Fy;Fz];
moment = [Mx;My;Mz];

% rotate to ground from contact point
fcp_g = W * force;
mcp_g = W * moment;

% Translate moments from tire_contact to wheel_center
fwc = fcp_g;
mwc = mcp_g + vec3_cross(wheel_center - tire_contact, fcp_g);

% Finally rotate from ground to wheel frame
wheel_center_force  = T * fwc;
wheel_center_moment = T * mwc;
end

function [rho, tire_contact, road_normal, road_grip, err_code] = ...
    road_model(jobflg, time, T, wheel_center, R_0, tire_width)

road_normal = [0;0;1];
road_origin = [0;0;0];
road_grip = 1.0;

% the contact point is the distance between the wheel_center and wheel_center along tire Z.
tire_Z = T(:,3);
loaded_radius = vec3_dot(wheel_center - road_origin, road_normal) / vec3_dot(tire_Z, road_normal);
tire_contact = wheel_center - loaded_radius * tire_Z; % vector from center to contact

rho = R_0 - loaded_radius;
rho = max(rho, 0.0);

err_code = 0; % no error
end

% Compute wheel center and slip velocities, Tydex-W orientation,
% slips kappa, alpha and inclination angle gamma.
function [Vc, Vs, W, kappa, alpha, gamma] = tire_kinematics(...
    frame_omega, T, wheel_center_vel, wheel_omega, rolling_radius, road_normal)

wheel_axis = T(:,2);

% W = transform from Tydex W tire frame to inertial frame
% W(:,3) == road surface normal, pointing upwards
%
W = zeros(3,3);
% Z axis: road surface normal, pointing upwards.
W(:,3) = road_normal;
% X axis: intersection of wheel plane with ground plane
W(:,1) = vec3_normalize(vec3_cross(wheel_axis, road_normal));
% Y axis: projection of wheel axis onto ground plane
W(:,2) = vec3_cross(road_normal, W(:,1));

% W' * T = rotate_x(gamma) = [1 0 0;
%                             0 cos(gamma) -sin(gamma);
%                             0 sin(gamma) cos(gamma)]:
% so gamma can be extracted from elements (3,2:3) of W' * T.
sin_gamma = vec3_dot(T(:,2), W(3,:));
cos_gamma = vec3_dot(T(:,3), W(3,:));
gamma = atan(sin_gamma / cos_gamma);

Vc = W' * wheel_center_vel;  % Tydex-W coordinates from inertial/road/ground coordinates

sign_Vcx = SGN(Vc(1));

Vcx_prime = Vc(1) + sign_Vcx*1e-6; % avoid division by zero by adding something small

% compute slip velocities
Vs = [-(Vc(1) - rolling_radius * wheel_omega); Vc(2)];

% compute standard slip quantities
kappa = Vs(1) / Vcx_prime;
tan_alpha = Vs(2) / Vcx_prime;

% Mirror alpha if traveling backwards
tan_alpha = tan_alpha * sign_Vcx;

alpha = atan(tan_alpha);
end

% "load" a tire.
% 1 = rear_mc - dunlop_F_100_90_19_D401.tir
% 2 = front_mc - dunlop_F_160_70_17_K591.tir
% 3 = tg2 - tg2.tir
function tire = load_tire(tire_id)

mm_to_m = 1/1000;

rear_mc.FILENAME = 'dunlop_R_160_70_17_K591.tir';
rear_mc.USE_MODE             = 2;
rear_mc.UNLOADED_RADIUS      = 320 * mm_to_m;
rear_mc.WIDTH                = 160.0 * mm_to_m;
rear_mc.VERTICAL_STIFFNESS   = 212.0 / mm_to_m; % N/mm -> N/m
rear_mc.VERTICAL_DAMPING     = 2.0 / mm_to_m; % N/(mm/s) -> N/(m/s)
rear_mc.ROLLING_RESISTANCE   = 0.02;
rear_mc.SHAPE = [
    1             0.0
    0.998312833   0.2
    0.993234135   0.4
    0.984711434   0.6
    0.972654199   0.8
    0.956928837   1.0
    0.680645161   0.711280855
    0.680645161   0.0];

% model parameters
rear_mc.CALPHA               = 1.70; % called CMX in .TIR file
rear_mc.CGAMMA               = 1.00;
rear_mc.UMIN                 = 1.40;
rear_mc.UMAX                 = 1.40;
rear_mc.RELAXATION_LENGTH  	= 100.0 * mm_to_m;
rear_mc.ALPHA_CRITICAL       = 6.0; % peak alpha
rear_mc.CURVATURE_FACTOR_ANGLE = 1.70;
rear_mc.SCALE_FACTOR_LATERAL = 0.715;
rear_mc.RATED_LOAD           = 1630;
rear_mc.SCALE_FACTOR_DIM     = -8E-5;
rear_mc.SLIP_RATIO_CRITICAL 	= 20.0;
rear_mc.CURVATURE_FACTOR_RATIO = 5.5;
rear_mc.PNEUM_TRAILING_SCALING = 0.5;
rear_mc.PNEUMATIC_LEAD_CAMBER  = 25.0;
rear_mc.LIMIT_CAMBER_ONSET_FRIC = 0.80;

front_mc.FILENAME = 'dunlop_F_100_90_19_D401.tir';
front_mc.USE_MODE            = 2;
front_mc.UNLOADED_RADIUS     = 341   * mm_to_m;
front_mc.WIDTH               = 100.0 * mm_to_m;
front_mc.VERTICAL_STIFFNESS  = 146.0 / mm_to_m; % N/mm -> N/m
front_mc.VERTICAL_DAMPING    = 0.2   / mm_to_m; % N/(mm/s) -> N/(m/s)
front_mc.ROLLING_RESISTANCE  = 0.02;
front_mc.CALPHA              = 1.70; % called CMX in .TIR file
front_mc.CGAMMA              = 1.00;
front_mc.UMIN                = 1.20; % UMIN > UMAX in .TIR file!
front_mc.UMAX                = 1.40;
front_mc.RELAXATION_LENGTH   = 100.0 * mm_to_m;
front_mc.ALPHA_CRITICAL      = 10.0; % peak alpha
front_mc.CURVATURE_FACTOR_ANGLE = 2.70;
front_mc.SCALE_FACTOR_LATERAL = 1.6417;
front_mc.RATED_LOAD          = 1662;
front_mc.SCALE_FACTOR_DIM    = -1.0E-4;
front_mc.SLIP_RATIO_CRITICAL = 20.0;
front_mc.CURVATURE_FACTOR_RATIO = 5.5;
front_mc.PNEUM_TRAILING_SCALING = 1.25;
front_mc.PNEUMATIC_LEAD_CAMBER 	= 20.0;
front_mc.LIMIT_CAMBER_ONSET_FRIC = 0.80;
front_mc.SHAPE = [
    1             0.0
    0.998312833   0.2
    0.993234135   0.4
    0.984711434   0.6
    0.972654199   0.8
    0.956928837   1.0
    0.680645161   0.711280855
    0.680645161   0.0 ];

tg2.FILENAME = 'tg2.tir';
tg2.USE_MODE            = 2;
tg2.UNLOADED_RADIUS     = 310 * mm_to_m;
tg2.WIDTH               = 180.0 * mm_to_m;
tg2.VERTICAL_STIFFNESS  = 230.0 / mm_to_m; % N/mm -> N/m
tg2.VERTICAL_DAMPING    = 2.0 / mm_to_m; % N/(mm/s) -> N/(m/s)
tg2.ROLLING_RESISTANCE  = 0.055;
tg2.CALPHA              = 2.20; % Called CMX in .TIR file
tg2.CGAMMA              = 1.40;
tg2.UMIN                = 1.30;
tg2.UMAX                = 1.50;
tg2.RELAXATION_LENGTH  	= 0.0 * mm_to_m; % no relaxation
tg2.ALPHA_CRITICAL      = 9.0;
tg2.CURVATURE_FACTOR_ANGLE = 1.85;
tg2.SCALE_FACTOR_LATERAL = 1.10;
tg2.RATED_LOAD          = 1000;
tg2.SCALE_FACTOR_DIM    = -1.5E-4;
tg2.SLIP_RATIO_CRITICAL = 20.0;
tg2.CURVATURE_FACTOR_RATIO = 5.5;
tg2.PNEUM_TRAILING_SCALING = 1.25;
tg2.PNEUMATIC_LEAD_CAMBER  = 25.0;
tg2.LIMIT_CAMBER_ONSET_FRIC = 0.30;
tg2.SHAPE = [
    1             0.0
    0.998312833   0.2
    0.993234135   0.4
    0.984711434   0.6
    0.972654199   0.8
    0.956928837   1.0
    0.680645161   0.711280855
    0.680645161   0.0];

tires = {front_mc, rear_mc, tg2};

tire = tires{tire_id};
end

% Compute loaded radius
% The tire is a linear spring, so this is easy.
function loaded_radius = loaded_radius_from_load(tire, Fz)
CFz = tire.VERTICAL_STIFFNESS;
unloaded_radius = tire.UNLOADED_RADIUS;
rho = Fz / CFz;

rho = max(rho, 0);  % loaded_radius always <= unloaded_radius.

loaded_radius = unloaded_radius - rho;
end

function test()

tire = load_tire(3);

Fz = 1500;
Vx = 100/3.6;
kappa0 = 0;
alpha0 = 0*pi/180;
curvature = 0;

% number of steps for continuous plots
steps = 150;

% number of plots in fugure
npv = 2;
nph = 2;

close all;
figure;

f = zeros(3,steps);
m = zeros(3,steps);

% Fy,Mz vs. alpha at gamma = 0 and gamma = 4
subplot(npv,nph,1);
hold on; grid on;

gamma = 0*pi/180;
alpha = linspace(-20,20,steps)*pi/180;
for i=1:steps
    [~,~,f(:,i),m(:,i)] = steady_tyr501(tire, Fz, Vx, kappa0, alpha(i), gamma, curvature);
end

plot(alpha*180/pi, f(2,:), 'b');
plot(alpha*180/pi, m(3,:), 'r');

gamma = 4*pi/180;
alpha = linspace(-20,20,steps)*pi/180;
for i=1:steps
    [~,~,f(:,i),m(:,i)] = steady_tyr501(tire, Fz, Vx, kappa0, alpha(i), gamma, curvature);
end

plot(alpha*180/pi, f(2,:), 'c');
plot(alpha*180/pi, m(3,:), 'm');
title('Fy,Mz vs. alpha at gamma=0 and 4');

% Fx vs kappa at alpha=gamma=0
subplot(npv,nph,2);
hold on; grid on;

gamma = 0*pi/180;
alpha = 0;
kappa = linspace(-1,1,steps);
for i=1:steps
    [~,~,f(:,i),m(:,i)] = steady_tyr501(tire, Fz, Vx, kappa(i), alpha, gamma, curvature);
end

plot(kappa, f(1,:), 'b');
title('Fx vs kappa');

% Fy vs. camber at alpha=kappa=0
subplot(npv,nph,3);
hold on; grid on;
alpha = 0;
kappa = 0;
gamma = linspace(-60,60,steps)*pi/180;
for i=1:steps
    [~,~,f(:,i),m(:,i)] = steady_tyr501(tire, Fz, Vx, kappa, alpha, gamma(i), curvature);
end

plot(gamma*180/pi, f(2,:), 'b');
title('Fy vs. camber at kappa=alpha=0');

subplot(npv,nph,4);
hold on; grid on;
axis equal;

% grip ellipse plot
nz = 8;
z_alpha = linspace(-16, 16, nz)*pi/180;
kappa = linspace(-1, 1, steps);
gamma = 0*pi/180;

fx = zeros(nz, steps);
fy = zeros(nz, steps);
mz = zeros(nz, steps);

colors = 'rgbgkmcrgbgkmcrgbgkmcrgbgkmc';

for iz = 1:nz
    for i=1:steps
        [~,~,f,m] = steady_tyr501(tire, Fz, Vx, kappa(i), z_alpha(iz), gamma, curvature);
        fx(iz,i) = f(1);
        fy(iz,i) = f(2);
        mz(iz,i) = m(3);
    end
    plot(fx(iz,:),fy(iz,:), [colors(iz) '.-']);
end

z_kappa = linspace(-.3, .3, nz);
alpha = linspace(-20, 20, steps)*pi/180;
for iz = 1:nz
    for i=1:steps
        [~,~,f,m] = steady_tyr501(tire, Fz, Vx, z_kappa(iz), alpha(i), gamma, curvature);
        fx(iz,i) = f(1);
        fy(iz,i) = f(2);
        mz(iz,i) = m(3);
    end
    plot(fx(iz,:),fy(iz,:), [colors(iz) '.-']);
end
title('grip ellipse');

end

function rolling_radius = compute_rolling_radius(tire, unclamped_nFz, R_omega, CFz)
Fz0 = tire.FNOMIN;
ndz = (Fz0/CFz);

rr_nFz = max(0.0, unclamped_nFz); % No Fz effect if wheel is off the ground
rolling_radius = R_omega - ndz * (tire.DREFF * atan(tire.BREFF*rr_nFz) + tire.FREFF*rr_nFz); % [Eqn (7) Paper]
end

% compute wheel orientation in inertial space from the wheel
% spin axis vector and the road normal.
function T = compute_wheel_orientation(wheel_axis, road_normal)
T = zeros(3,3);
T(:,1) = vec3_normalize(vec3_cross(wheel_axis, road_normal));
T(:,2) = wheel_axis;
T(:,3) = vec3_cross(T(:,1), T(:,2));
end

% fwc,mwc are forces and moments at hub, in the inertial/ground/road frame.
% fcp,mcp are forces at contact point in Tydex-W wheel frame.
function [fwc, mwc, fcp, mcp] = steady_tyr501(tire, Fz, Vx, kappa, alpha, gamma, curvature)

iswtch = 1; % 0 = static, 1 = dynamic
jobflg = 0; %normal;

Vcx = Vx;

loaded_radius = loaded_radius_from_load(tire, Fz);

unloaded_radius = tire.UNLOADED_RADIUS;

% In lieu of a real rolling radius model, we will
% use the constant unloaded radius.
rolling_radius = unloaded_radius;

road_normal = [0;0;1];
wheel_axis = [sin(alpha)*cos(gamma); cos(alpha)*cos(gamma); sin(gamma)];

T = compute_wheel_orientation(wheel_axis, road_normal);

% Compute contact patch location
radial = loaded_radius * T(:,3);
tire_contact = [radial(1:2);0];

wheel_center     = [0;0;radial(3)];
wheel_center_vel = [Vcx/cos(alpha);0;0];   % wheel center velocity

wheel_omega = (kappa + 1) * Vcx / rolling_radius;

frame_omega = [0;0;curvature*Vcx];

time = 1.0; % set time past startup threshold so we get some forces

wheel_angle = -1;

% these integration states are just the slips.
% Initialize so they aren't lagged
states = [kappa;alpha];
[fwc, mwc, fcp, mcp, states_dot] = ...
    tyr501Harty_2020( tire, iswtch, jobflg, ...
    time, ...,
    frame_omega, ...
    wheel_center, T, wheel_center_vel, ...
    wheel_angle, wheel_omega, ...
    states);
end

% Handy functions for smooth transition between s = 0 and s = 1.
% smoothstep is a cubic with first derivatives == 0 at ends, y = 0 at s = 0,
% and y = 1 at s = 1.  The quintic smootherstep is similar,
% but both first and second derivatives are 0 at ends.
% smoothstep is closest to (1-cos(pi*s))/2; smootherstep is more gradual at ends,
% steeper at s = 0.5.  Both are much more efficient than (1-cos(pi*s))/2.

function s = smoothstep(s)
s = (s*s*(3 - 2*s));
end

function s = smootherstep(s)
s = s*s*s*(s*(s * 6 - 15) + 10);
end

% step() interpolates smoothly between y0 and y1
% depending on the value of x.
% y0 is returned if x <= x0, y1 if x >= x1.
% In between step smoothly interpolated between
% y0 and y1, with derivative 0 at x=x0 and x=x1.
% Similar to 1-cos((x-x0)/(x1-x0)*pi).
function y = step(x, x0,y0, x1,y1)
if (x < x0)
    y = y0;
elseif (x > x1)
    y = y1;
else
    s = (x - x0) / (x1 - x0);
    s = smoothstep(s);
    y = y0 + (y1 - y0) * s;
end
end

% Returns -1 if less than zero, 1 otherwise
% (unlike sign() function which returns 0 if 0)
function s = SGN(x)
s = 1;
if (x < 0), s = -1; end
end

% Handy 3-vector functions
function n = vec3_norm(v)
n = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
end

function v = vec3_normalize(v)
n = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
% norm assumed to be > 0
%assert(n > 0);
v = v / n;
end

function d = vec3_dot(a, b)
d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
end

function c = vec3_cross(a, b)
c = [a(2)*b(3) - a(3)*b(2);
    a(3)*b(1) - a(1)*b(3);
    a(1)*b(2) - a(2)*b(1)];
end
