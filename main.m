%% ========================================================================
%% PROJECT : Interplanetary 3D Low-Thrust Transfer Optimization (ICLOCS2 + IPOPT)
%  AUTHOR  : Gabriele Borghini - November 2026
%  Based on ICLOCS2 problem template (Copyright (c) 2018 Imperial College London, MIT License)
% =========================================================================
% COPYRIGHT AND USAGE RIGHTS:
%  The source code, problem formulation, and custom scripts contained herein 
%  are the intellectual property of the author. They may not be copied, 
%  modified, redistributed, or used for commercial purposes without the 
%  express written consent of the author.
%
% ACADEMIC CONTEXT & THIRD-PARTY TOOLS:
%  Developed solely for academic research, advanced astrodynamics studies, 
%  and educational purposes at Alma Mater Studiorum - University of Bologna.
%  This implementation interfaces with ICLOCS2 (Imperial College London 
%  Optimal Control Software) and IPOPT, which remain the property of their 
%  respective authors and are subject to their own licensing terms.
%
% DISCLAIMER OF WARRANTY:
%  This software is provided "as is", without warranty of any kind, express
%  or implied, including but not limited to the warranties of accuracy, 
%  numerical stability, or fitness for a particular purpose. In no event 
%  shall the author be liable for any claim or damages arising from the use 
%  of this software.
% -------------------------------------------------------------------------
% DESCRIPTION:
%  This code is an extension and adaptation of the direct collocation framework 
%  originally developed for launcher ascent trajectory optimization (using ICLOCS2 
%  and IPOPT). The formulation is refactored to solve 3D heliocentric low-thrust 
%  interplanetary transfers with a prescribed arrival epoch, free departure time 
%  (free TOF), and Cartesian steering controls [ux, uy, uz, tau]. 
%
%  The script automatically extracts Earth and target ephemeris data via vector 
%  tables, evaluates continuous boundary conditions through smooth polynomial 
%  approximations, and employs a homotopy continuation method (energy-optimal 
%  to fuel-optimal) to ensure robust convergence to the optimal bang-bang-like 
%  trajectory.
% 
%  Final validation and open-loop verification of the optimal 
%  solution are performed via MATLAB's ode45 integrator with strict tolerances 
%  (RelTol/AbsTol less or equal to 1e-8).
% =========================================================================
%  Copyright (c) 2026 Gabriele Borghini. All Rights Reserved.
% =========================================================================

clc; clear; close all; format compact;

if exist('ephemeris','dir'),       addpath(genpath('ephemeris'));
elseif exist('../ephemeris','dir'), addpath(genpath('../ephemeris')); end
addpath(genpath('D:\OptimalControl\ICLOCS-master\src'))
addpath(genpath('D:\OptimalControl\ICLOCS-master\usr'))
addpath(genpath('D:\OptimalControl\adigator-master'));
delete('Func_*_low_thrust_dynamics.m');

%% --- HELIOCENTRIC CANONICAL UNITS ---
AU     = 1.495978707e8;                 % [km] Length unit (LU)
mu_sun = 1.32712440018e11;              % [km^3/s^2]
VU     = sqrt(mu_sun / AU);             % [km/s] Velocity unit (~29.7847 km/s)
TU     = AU / VU;                       % [s] Time unit (~5.0226e6 s = 58.1324 days)
g0     = 9.80665e-3;                    % [km/s^2]

%% --- VEHICLE DATA (e.g., Hall / Ion Thruster) ---
m0_kg    = 1000.0;                      % [kg] Initial mass (MU)
m_dry_kg = 600.0;                       % [kg] Minimum dry mass
Thrust_N = 0.35;                        % [N] Maximum thrust (350 mN)
Isp_s    = 3000.0;                      % [s] Specific impulse

% Conversion to canonical dimensionless units
Thrust_km = Thrust_N * 1e-3;                  % [kg * km / s^2]
mdot_kgs  = Thrust_km / (Isp_s * g0);         % [kg/s]

Tmax_nd  = (Thrust_km / m0_kg) * (TU^2 / AU); % Initial acceleration in [AU/TU^2]
mdot_nd  = (mdot_kgs / m0_kg) * TU;           % Mass flow rate in [MU/TU]
m_dry_nd = m_dry_kg / m0_kg;                  % Dimensionless dry mass [-]

%% --- BOUNDARY CONDITIONS EXTRACTION (FIXED ARRIVAL TIME, FREE DEPARTURE) ---
file_earth  = 'ephem_earth.txt';        % Your Earth ephemeris file
file_target = 'ephem_target.txt';       % Your Target ephemeris file

tf_epoch        = '2029-03-25 00:00:00'; % FIXED Target Arrival Date!
TOF_guess_days  = 460.0;                 % [days] Initial TOF guess
TOF_window_days = [360.0, 480.0];        % [days] Window [TOF_min, TOF_max] -> optimizes departure date!
n_rev           = [];                    % Automatic computation of revolutions
units_flag      = 'km_s';

[x0_guess_nd, xf_nd, tf_bounds_tu, ephem_info] = get_boundary_conditions( ...
    file_earth, file_target, tf_epoch, TOF_guess_days, TOF_window_days, n_rev, units_flag);

rf = xf_nd(1);  thetaf = xf_nd(2);  phif = xf_nd(3);  vf = xf_nd(4);

fprintf('Departure Guess (TOF %.0f d) -> r0 = %.4f AU, lon0 = %.2f deg, lat0 = %.2f deg (TOF Window: [%.0f, %.0f] d)\n', ...
        TOF_guess_days, x0_guess_nd(1), rad2deg(x0_guess_nd(2)), rad2deg(x0_guess_nd(3)), TOF_window_days(1), TOF_window_days(2));
fprintf('Fixed Arrival (JD %.2f)   -> rf = %.4f AU, lonf = %.2f deg, latf = %.2f deg\n', ...
        ephem_info.JD_tf, rf, rad2deg(thetaf), rad2deg(phif));

%% --- PROBLEM (STEP 1: ENERGY-OPTIMAL -> STEP 2: FUEL-OPTIMAL) ---
N_nodes_step1 = 150;  % Grid for Step 1 (energy)
N_nodes_hom   = 300;  % Dense grid for energy -> fuel homotopy (exact match with ode45)

% Residual checks both Earth departure matching at t=0 and Target arrival matching at t=tf
res_fun = @(sol) max([max(abs(sol.X(1,1:6).'   - ppval(ephem_info.pp_earth, sol.tf))); ...
                      max(abs(sol.X(end,1:6).' - xf_nd(1:6).'))]);
res_tol = 5e-3;

fprintf('n_rev = %d (estimate %d), margin = %.2f\n', ephem_info.n_rev, ephem_info.n_rev_est, ephem_info.dtheta_margin);

fprintf('\n=== STEP 1: Energy-Optimal Resolution with free departure (free TOF) ===\n');
[problem_e, guess_e] = low_thrust(x0_guess_nd, xf_nd, tf_bounds_tu, ephem_info.pp_earth, Tmax_nd, mdot_nd, m_dry_nd, 'energy');
options_e = problem_e.settings(N_nodes_step1);
options_e.meshstrategy = 'fixed';
[sol_e, ~] = solveMyProblem(problem_e, guess_e, options_e);

r1 = res_fun(sol_e);
fprintf('Step 1: boundary residual = %.2e (tol %.0e)\n', r1, res_tol);
if r1 > res_tol
    error(['Step 1 infeasible (residual %.2e). Cannot proceed: widen TOF_window_days, ' ...
           'check n_rev and maximum thrust.'], r1);
end

fprintf('\n=== STEP 2: Homotopy energy -> fuel, free departure (free TOF) ===\n');
eps_list = [1e-2, 5e-4, 2e-5];
sol_prev = sol_e;
solution = sol_e;
problem  = problem_e;

for k = 1:numel(eps_list)
    fprintf('\n--- Homotopy: eps = %.1e ---\n', eps_list(k));
    [problem_k, guess_k] = low_thrust(x0_guess_nd, xf_nd, tf_bounds_tu, ephem_info.pp_earth, Tmax_nd, mdot_nd, m_dry_nd, eps_list(k));
    guess_k.tf          = sol_prev.tf;
    guess_k.time        = sol_prev.T;
    guess_k.states      = sol_prev.X;
    guess_k.inputs      = sol_prev.U;
    if isfield(sol_prev, 'multipliers')
        guess_k.multipliers = sol_prev.multipliers;
    end
    options_k = problem_k.settings(N_nodes_hom);
    options_k.meshstrategy = 'fixed';
    if k > 1
        options_k.start = 'Warm';
        options_k.ipopt.mu_init = 1e-4;
    end
    [sol_k, ~] = solveMyProblem(problem_k, guess_k, options_k);

    rk = res_fun(sol_k);
    fprintf('eps = %.1e: boundary residual = %.2e | m_f = %.2f kg\n', eps_list(k), rk, sol_k.X(end,7)*m0_kg);
    if rk > res_tol
        warning('Homotopy interrupted at eps = %.1e (residual %.2e). Keeping previous solution.', eps_list(k), rk);
        break
    end
    sol_prev = sol_k;  solution = sol_k;  problem = problem_k;
end

%% --- RESULTS & ODE45 VERIFICATION ---
vd   = problem.data;
Tn   = solution.T(:);  
Un   = solution.U;
uRaw = @(t) interp1(Tn, Un, t, 'linear', 'extrap');
uFun = @(t) [uRaw(t)*[eye(3); zeros(1,3)] ./ max(sqrt(sum((uRaw(t)*[eye(3); zeros(1,3)]).^2, 2)), 1e-12), ...
             max(0, min(1, uRaw(t)*[0;0;0;1]))];
rhs  = @(t,x) low_thrust_dynamics(x.', uFun(t), [], t, vd).';

odeopt = odeset('RelTol', 1e-10, 'AbsTol', 1e-10*ones(1,7));
[tv, xv] = ode45(rhs, [Tn(1) Tn(end)], solution.X(1,:).', odeopt);
uv = uFun(tv);

%% --- DATA REFINEMENT ---
t_days_opt   = solution.T * TU / 86400;
t_days_sim   = tv * TU / 86400;
TOF_opt_days = solution.tf * TU / 86400;
JD_dep_opt   = ephem_info.JD_tf - TOF_opt_days;

% Evaluate exact Earth state at the optimal departure epoch found by IPOPT
x0_earth_opt = ppval(ephem_info.pp_earth, solution.tf);
r0     = x0_earth_opt(1);
theta0 = x0_earth_opt(2);
phi0   = x0_earth_opt(3);
v0     = x0_earth_opt(4);

mf_opt_kg  = solution.X(end,7) * m0_kg;
mf_sim_kg  = xv(end,7) * m0_kg;
mp_used_kg = m0_kg - mf_opt_kg;
dV_tot_kms = Isp_s * g0 * log(m0_kg / mf_opt_kg);

%% --- FINAL REPORT ---
fprintf('\n=== LOW-THRUST TRANSFER RESULTS (FIXED ARRIVAL, FREE DEPARTURE) ===\n\n');
fprintf('OPTIMAL time of flight (TOF): %.2f days (Guess: %.1f d | Window: [%.0f, %.0f] d)\n', ...
        TOF_opt_days, TOF_guess_days, TOF_window_days(1), TOF_window_days(2));
fprintf('OPTIMAL departure date      : %s (JD %.2f)\n', ...
        datestr(datetime(JD_dep_opt, 'ConvertFrom', 'juliandate')), JD_dep_opt);
fprintf('FIXED arrival date          : %s (JD %.2f)\n', ...
        datestr(datetime(ephem_info.JD_tf, 'ConvertFrom', 'juliandate')), ephem_info.JD_tf);
fprintf('Initial radius (opt/Earth)  : %.4f AU / %.4f AU\n', solution.X(1,1), r0);
fprintf('Final radius (opt / targ)   : %.4f AU / %.4f AU (ode45: %.4f AU)\n', solution.X(end,1), rf, xv(end,1));
fprintf('Final longitude             : %.2f deg (target: %.2f deg | ode45: %.2f deg)\n', rad2deg(solution.X(end,2)), rad2deg(thetaf), rad2deg(xv(end,2)));
fprintf('Final latitude              : %.2f deg (target: %.2f deg | ode45: %.2f deg)\n', rad2deg(solution.X(end,3)), rad2deg(phif), rad2deg(xv(end,3)));
fprintf('Final velocity              : %.4f km/s (target: %.4f km/s | ode45: %.4f km/s)\n', solution.X(end,4)*VU, vf*VU, xv(end,4)*VU);
fprintf('Final mass (IPOPT)          : %.2f kg (Propellant consumed: %.2f kg)\n', mf_opt_kg, mp_used_kg);
fprintf('Final mass (ode45)          : %.2f kg\n', mf_sim_kg);
fprintf('Equivalent total Delta-V    : %.3f km/s\n', dV_tot_kms);

%% --- STATE EXTRACTION AND KINEMATIC PLOTS ---
r_opt     = solution.X(:,1);
theta_opt = solution.X(:,2);
phi_opt   = solution.X(:,3);
v_opt     = solution.X(:,4) * VU;         % [km/s]
gamma_opt = solution.X(:,5);
chi_opt   = solution.X(:,6);
m_opt_kg  = solution.X(:,7) * m0_kg;      % [kg]

ux_opt       = solution.U(:,1);
uy_opt       = solution.U(:,2);
uz_opt       = solution.U(:,3);
throttle_opt = solution.U(:,4);

% Reconstruct equivalent pointing angles alpha (pitch) and beta - chi (in-plane offset)
alpha_opt    = asin(max(-1.0, min(1.0, ux_opt)));
uh_opt       = uy_opt .* cos(chi_opt) + uz_opt .* sin(chi_opt);
uc_opt       = uz_opt .* cos(chi_opt) - uy_opt .* sin(chi_opt);
beta_rel_opt = atan2(uc_opt, uh_opt);

figure('Name', 'Low-Thrust Kinematics & Controls', 'Color', 'w');

% 1. Heliocentric Distance
subplot(2, 2, 1); hold on;
plot(t_days_opt, r_opt, 'r', 'LineWidth', 2);
plot(t_days_sim, xv(:,1), 'k--', 'LineWidth', 1.5);
grid on; xlabel('Time since Departure [days]'); ylabel('Heliocentric Distance r [AU]');
title('Heliocentric Radius vs Time');
legend('IPOPT', 'ode45', 'Location', 'best'); theme('light');

% 2. Throttling Profile and Mass
subplot(2, 2, 2);
yyaxis left
area(t_days_opt, throttle_opt, 'FaceColor', [0.2 0.6 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
ylabel('Throttle \tau [-]'); ylim([-0.05 1.1]);
yyaxis right
plot(t_days_opt, m_opt_kg, 'LineWidth', 2, 'Color', [0.85 0.33 0.10]);
ylabel('Spacecraft Mass [kg]');
grid on; xlabel('Time since Departure [days]');
title('Thrust Throttle Profile & Spacecraft Mass'); theme('light');

% 3. Pointing Angles (Reconstructed from ux, uy, uz)
subplot(2, 2, 3); hold on;
plot(t_days_opt, rad2deg(alpha_opt), 'm', 'LineWidth', 1.8);
plot(t_days_opt, rad2deg(beta_rel_opt), 'c-', 'LineWidth', 1.8);
grid on; xlabel('Time since Departure [days]'); ylabel('Control Angles [deg]');
legend('Out-of-plane / Pitch \alpha', 'In-plane offset (\beta - \chi)', 'Location', 'best');
title('Thrust Pointing Angles vs Time'); theme('light');

% 4. Velocity Components (Radial, Transverse, Normal)
v_rad   = v_opt .* sin(gamma_opt);
v_trans = v_opt .* cos(gamma_opt) .* sin(chi_opt);
v_norm  = v_opt .* cos(gamma_opt) .* cos(chi_opt);

subplot(2, 2, 4); hold on;
plot(t_days_opt, v_rad, 'b', 'LineWidth', 1.8);
plot(t_days_opt, v_trans, 'g', 'LineWidth', 1.8);
plot(t_days_opt, v_norm, 'r', 'LineWidth', 1.8);
grid on; xlabel('Time since Departure [days]'); ylabel('Velocity [km/s]');
legend('v_r (Radial)', 'v_\theta (Transverse)', 'v_\phi (Normal)', 'Location', 'best');
title('Heliocentric Velocity Components vs Time'); theme('light');

% Additional Plot: Cartesian Controls ux, uy, uz & Norm (as in ascent main.m)
figure('Name', 'Cartesian Thrust Direction Controls vs Time', 'Color', 'w'); hold on;
u_dir_norm = sqrt(ux_opt.^2 + uy_opt.^2 + uz_opt.^2);
plot(t_days_opt, ux_opt, 'b', 'LineWidth', 2);
plot(t_days_opt, uy_opt, 'r', 'LineWidth', 2);
plot(t_days_opt, uz_opt, 'g', 'LineWidth', 2);
plot(t_days_opt, u_dir_norm, 'k:', 'LineWidth', 1.5);
grid on; xlabel('Time since Departure [days]'); ylabel('Control components [-]');
legend('u_x (Up / Radial)', 'u_y (North)', 'u_z (East)', '||u_{dir}||', 'Location', 'best');
title('Thrust Direction Controls (u_x, u_y, u_z) vs Time'); theme('light');

%% --- COSTATES (CARTESIAN), SWITCHING FUNCTION & HAMILTONIAN ---
if isfield(solution, 'multipliers') && isfield(solution.multipliers, 'lambda_1toN')
    lambda_sph = solution.multipliers.lambda_1toN;
    N_opt = length(t_days_opt);
    
    lambda_r_cart = zeros(N_opt, 3);
    lambda_v_cart = zeros(N_opt, 3);
    H_tot         = zeros(N_opt, 1);
    SF            = zeros(N_opt, 1); % Switching Function
    
    for k = 1:N_opt
        rk   = r_opt(k);
        thk  = theta_opt(k);
        phik = phi_opt(k);
        vk   = max(solution.X(k,4), 1e-6);
        gamk = gamma_opt(k);
        chik = chi_opt(k);
        mk   = solution.X(k,7);
        
        lr   = lambda_sph(k, 1);
        lth  = lambda_sph(k, 2);
        lphi = lambda_sph(k, 3);
        lv   = lambda_sph(k, 4);
        lgam = lambda_sph(k, 5);
        lchi = lambda_sph(k, 6);
        lm   = lambda_sph(k, 7);
        
        cos_gam = cos(gamk);
        if abs(cos_gam) < 1e-6; cos_gam = 1e-6; end
        
        % --- A) Velocity costates ---
        lv_U = lv * sin(gamk) + (lgam / vk) * cos(gamk);
        lv_N = (lv * cos(gamk) - (lgam / vk) * sin(gamk)) * cos(chik) - (lchi / (vk * cos_gam)) * sin(chik);
        lv_E = (lv * cos(gamk) - (lgam / vk) * sin(gamk)) * sin(chik) + (lchi / (vk * cos_gam)) * cos(chik);
        
        R_loc2heci = [cos(phik)*cos(thk), -sin(phik)*cos(thk), -sin(thk);
                      cos(phik)*sin(thk), -sin(phik)*sin(thk),  cos(thk);
                      sin(phik),           cos(phik),           0];
                  
        lambda_v_cart(k, :) = (R_loc2heci * [lv_U; lv_N; lv_E]).';
        
        % --- B) Position costates ---
        v_U = vk * sin(gamk);
        v_N = vk * cos(gamk) * cos(chik);
        v_E = vk * cos(gamk) * sin(chik);
        cos_phi = max(abs(cos(phik)), 1e-6) * sign(cos(phik) + 1e-12);
        
        lr_U = lr;
        lr_N = (lphi - (lv_U * v_N - lv_N * v_U)) / rk;
        lr_E = (lth / cos_phi - lv_U * v_E + lv_E * v_U + (lv_N * v_E - lv_E * v_N) * tan(phik)) / rk;
        
        lambda_r_cart(k, :) = (R_loc2heci * [lr_U; lr_N; lr_E]).';
        
        % --- C) Primer Vector norm and Switching Function ---
        p_norm = norm(lambda_v_cart(k, :));
        SF(k)  = (p_norm * Tmax_nd / mk) - abs(lm) * mdot_nd;
        
        % --- D) Complete Hamiltonian (H = L + lambda^T * f) ---
        dx_k  = low_thrust_dynamics(solution.X(k,:), solution.U(k,:), [], solution.T(k), problem.data);
        eps_f = problem.data.mode_opt;
        if ischar(eps_f), eps_f = 0; end
        L_k   = 10.0 * eps_f * (solution.U(k,4)^2) + 1e-7 * sum(solution.U(k,1:3).^2);
        H_tot(k) = L_k + dot(lambda_sph(k,:), dx_k(:).');
    end
    
    figure('Name', 'Costates, Primer Vector & Hamiltonian', 'Color', 'w');
    subplot(3,1,1); hold on;
    plot(t_days_opt, lambda_v_cart(:,1), 'r', 'LineWidth', 1.5);
    plot(t_days_opt, lambda_v_cart(:,2), 'g', 'LineWidth', 1.5);
    plot(t_days_opt, lambda_v_cart(:,3), 'b', 'LineWidth', 1.5);
    plot(t_days_opt, vecnorm(lambda_v_cart, 2, 2), 'k--', 'LineWidth', 1.8);
    grid on; ylabel('\lambda_v [-]');
    legend('\lambda_{vx}', '\lambda_{vy}', '\lambda_{vz}', '||\lambda_v|| (Primer Mag)', 'Location', 'best');
    title('Cartesian Velocity Costates (Primer Vector)'); theme('light');
    
    subplot(3,1,2); hold on;
    plot(t_days_opt, SF, 'b', 'LineWidth', 1.8);
    yline(0, 'r--', 'LineWidth', 1.2);
    grid on; ylabel('Switching Func.');
    title('Switching Function vs Time'); theme('light');
    
    subplot(3,1,3); hold on;
    plot(t_days_opt, H_tot, 'k', 'LineWidth', 1.8);
    grid on; xlabel('Time since Departure [days]'); ylabel('Hamiltonian H');
    title('Hamiltonian vs Time (Autonomous System => H = const)'); theme('light');
end

%% --- 3D HELIOCENTRIC TRAJECTORY PLOT (HECI) ---
X_opt = r_opt .* cos(phi_opt) .* cos(theta_opt);
Y_opt = r_opt .* cos(phi_opt) .* sin(theta_opt);
Z_opt = r_opt .* sin(phi_opt);

% Direct transformation of local Cartesian thrust [ux, uy, uz] * throttle to HECI frame
u_loc_U = throttle_opt .* ux_opt;
u_loc_N = throttle_opt .* uy_opt;
u_loc_E = throttle_opt .* uz_opt;

Tx_heci = u_loc_U.*cos(phi_opt).*cos(theta_opt) - u_loc_N.*sin(phi_opt).*cos(theta_opt) - u_loc_E.*sin(theta_opt);
Ty_heci = u_loc_U.*cos(phi_opt).*sin(theta_opt) - u_loc_N.*sin(phi_opt).*sin(theta_opt) + u_loc_E.*cos(theta_opt);
Tz_heci = u_loc_U.*sin(phi_opt)                 + u_loc_N.*cos(phi_opt);

figure('Name', '3D Heliocentric Low-Thrust Trajectory', 'Color', 'w'); hold on;

th_ring = linspace(0, 2*pi, 200);
plot3(r0*cos(th_ring), r0*sin(th_ring), zeros(size(th_ring)), 'b:', 'LineWidth', 1.2);
plot3(rf*cos(th_ring), rf*sin(th_ring), zeros(size(th_ring)), 'r:', 'LineWidth', 1.2);

plot3(0, 0, 0, 'yo', 'MarkerFaceColor', [1 0.85 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
plot3(X_opt, Y_opt, Z_opt, 'k', 'LineWidth', 2.2);
plot3(X_opt(1), Y_opt(1), Z_opt(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot3(X_opt(end), Y_opt(end), Z_opt(end), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 10);

idx_q = 1:6:length(X_opt);
quiver3(X_opt(idx_q), Y_opt(idx_q), Z_opt(idx_q), ...
        Tx_heci(idx_q), Ty_heci(idx_q), Tz_heci(idx_q), 0.35, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);

grid on; axis equal;
xlabel('X_{HECI} [AU]'); ylabel('Y_{HECI} [AU]'); zlabel('Z_{HECI} [AU]');
title('3D Heliocentric Low-Thrust Transfer & Thrust Vectors');
legend('Earth Orbit (1 AU)', 'Target Orbit', 'Sun', 'Transfer Trajectory', ...
       'Optimal Earth Departure (t_0)', 'Fixed Target Arrival (t_f)', 'Thrust Vectors', 'Location', 'best');
view(35, 25); theme('light');