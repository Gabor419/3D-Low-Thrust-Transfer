function [problem, guess] = low_thrust(x0_guess_nd, xf_guess_nd, tf_bounds_tu, dtf_bounds_tu, ephem_info, Tmax_nd, mdot_nd, m_dry_nd, mode_opt)
%% ========================================================================
%% PROJECT : Interplanetary 3D Low-Thrust Transfer Optimization (ICLOCS2 + IPOPT)
%  AUTHOR  : Gabriele Borghini - November 2026
%  Based on ICLOCS2 problem template (Copyright (c) 2018 Imperial College London, MIT License)
% =========================================================================
%  Copyright (c) 2026 Gabriele Borghini. All Rights Reserved.
% 
% COPYRIGHT AND USAGE RIGHTS:
%  The source code, problem formulation, and custom scripts contained herein 
%  are the intellectual property of the author. Permission is hereby granted,
%  free of charge, to any person obtaining a copy of this software and 
%  associated documentation files, to deal in the Software without 
%  restriction, including without limitation the rights to use, copy, 
%  modify, merge, publish, distribute and sublicense copies of the Software,
%  and to permit persons to whom the Software is furnished to do so, subject 
%  to the following conditions:
%  The above copyright notice and this permission notice shall be included 
%  in all copies or substantial portions of the Software.
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
%  3D heliocentric low-thrust interplanetary transfer (Fixed range Arrival, Free Departure).
%   - fixed arrival range at target body; free TOF (free departure epoch)
%   - continuous departure matching via C^inf smooth Earth ephemeris polynomials
%   - Cartesian unit vector steering controls [ux, uy, uz] and throttling tau in [0,1]
%   - cost function: fuel-optimal maximization of final mass (-m_f) via homotopy
%
%  State:   x = [r, theta, phi, v, gamma, chi, m]
%  Control: u = [ux, uy, uz, tau]   (Up, North, East, Throttle)
%  Normalization Units: Canonical heliocentric units (AU, TU, VU, MU)
% =========================================================================

%% --- CONSTANTS ---
mu = 1.0;

%% --- MODEL ---
InternalDynamics = @low_thrust_dynamics;
SimDynamics      = @low_thrust_dynamics;

problem.analyticDeriv.gradCost          = [];
problem.analyticDeriv.hessianLagrangian = [];
problem.analyticDeriv.jacConst          = [];
problem.settings = @settings_low_thrust;

%% --- TIME & PARAMETERS (Clock [0, tf] where tf = TOF is FREE, p(1) = dtf is FREE) ---
problem.time.t0_min = 0;                  problem.time.t0_max = 0;                  guess.t0 = 0;
problem.time.tf_min = tf_bounds_tu.min;   problem.time.tf_max = tf_bounds_tu.max;   guess.tf = tf_bounds_tu.guess;

problem.parameters.pl = dtf_bounds_tu.min;
problem.parameters.pu = dtf_bounds_tu.max;
guess.parameters      = dtf_bounds_tu.guess;

%% --- STATES ---
theta_min = min(x0_guess_nd(2), xf_guess_nd(2)) - 2*pi;
theta_max = max(x0_guess_nd(2), xf_guess_nd(2)) + 2*pi;
ang_lim   = 1.2;

% Path bounds
problem.states.xl = [0.3,  theta_min, -ang_lim,  0.2,  -ang_lim,  -pi,   m_dry_nd];
problem.states.xu = [4.0,  theta_max,  ang_lim,  2.5,   ang_lim,   pi,   1.0];

% Initial state: x0(1:6) is FREE within bounds (matched to Earth via b_unscaled), m(0) = 1.0 FIXED
problem.states.x0  = [];
problem.states.x0l = [problem.states.xl(1:6), 1.0];
problem.states.x0u = [problem.states.xu(1:6), 1.0];

% Final state: xf(1:6) is FREE within bounds (matched to Target via b_unscaled), m(tf) is FREE in [m_dry_nd, 1.0]
problem.states.xfl = [problem.states.xl(1:6), m_dry_nd];
problem.states.xfu = [problem.states.xu(1:6), 1.0];

problem.states.xErrorTol_local    = [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6];
problem.states.xErrorTol_integral = [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6];
problem.states.xConstraintTol     = [1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5];

%% --- INPUT: u = [ux, uy, uz, throttle] ---
problem.inputs.N   = 0;
problem.inputs.ul  = [-1.1, -1.1, -1.1, 0.0];
problem.inputs.uu  = [ 1.1,  1.1,  1.1, 1.0];
problem.inputs.rl  = [-2.0, -2.0, -2.0, -inf];
problem.inputs.ru  = [ 2.0,  2.0,  2.0,  inf];
problem.inputs.u0l = [-1.1, -1.1, -1.1, 0.0];
problem.inputs.u0u = [ 1.1,  1.1,  1.1, 1.0];
problem.inputs.uConstraintTol = [1e-4, 1e-4, 1e-4, 1e-4];

%% --- PATH CONSTRAINTS AND BOUNDARY CONSTRAINTS ---
problem.setpoints.states = [];
problem.setpoints.inputs = [];

% Equality constraint (ng_eq = 1): ux^2 + uy^2 + uz^2 - 1 = 0
problem.constraints.ng_eq    = 1;
problem.constraints.gTol_eq  = 1e-5;
problem.constraints.gl       = [];
problem.constraints.gu       = [];
problem.constraints.gTol_neq = [];

% Boundary constraint b_unscaled (12 equations):
%  1..6:  x0(1:6) - x_Earth(dt0)  = 0,  with dt0 = p(1) - (tf - tf_guess)
%  7..12: xf(1:6) - x_Target(dtf) = 0,  with dtf = p(1)
problem.constraints.bl   = zeros(12, 1);
problem.constraints.bu   = zeros(12, 1);
problem.constraints.bTol = 1e-5 * ones(12, 1);

%% --- DATA ---
problem.data.mu          = mu;
problem.data.Tmax        = Tmax_nd;
problem.data.mdot_max    = mdot_nd;
problem.data.mode_opt    = mode_opt;
problem.data.pp_earth    = ephem_info.pp_earth;
problem.data.pp_target   = ephem_info.pp_target;
problem.data.P_earth     = ephem_info.poly_earth;
problem.data.P_target    = ephem_info.poly_target;
problem.data.dt0_mid     = ephem_info.dt0_mid;
problem.data.dt0_half    = ephem_info.dt0_half;
problem.data.dtf_mid     = ephem_info.dtf_mid;
problem.data.dtf_half    = ephem_info.dtf_half;
problem.data.tf_guess_tu = tf_bounds_tu.guess;

%% --- SMOOTH INITIAL GUESS (4 Controls: ux, uy, uz, throttle) ---
N_guess = 100;
t_vec   = linspace(0, tf_bounds_tu.guess, N_guess).';
s       = linspace(0, 1, N_guess).';
guess.time = t_vec;

r_g     = (1 - s)*x0_guess_nd(1) + s*xf_guess_nd(1);
theta_g = (1 - s)*x0_guess_nd(2) + s*xf_guess_nd(2);
phi_g   = (1 - s)*x0_guess_nd(3) + s*xf_guess_nd(3);
v_g     = (1 - s)*x0_guess_nd(4) + s*xf_guess_nd(4);
gam_g   = (1 - s)*x0_guess_nd(5) + s*xf_guess_nd(5);
chi_g   = (1 - s)*x0_guess_nd(6) + s*xf_guess_nd(6);
m_tf_g  = max(m_dry_nd + 0.1, 1.0 - 0.3 * mdot_nd * tf_bounds_tu.guess);
m_g     = (1 - s)*1.0 + s*m_tf_g;

guess.states = [r_g, theta_g, phi_g, v_g, gam_g, chi_g, m_g];

ux_g  = zeros(N_guess, 1);
uy_g  = cos(chi_g);
uz_g  = sin(chi_g);
thr_g = 0.5 * ones(N_guess, 1);

guess.inputs = [ux_g, uy_g, uz_g, thr_g];

%% --- HANDLES (Standard ICLOCS2) ---
problem.data.InternalDynamics   = InternalDynamics;
problem.data.functionfg         = @fg;
problem.data.plantmodel         = func2str(InternalDynamics);
problem.functions               = {@L, @E, @f, @g, @avrc, @b};
problem.sim.functions           = SimDynamics;
problem.sim.inputX              = [];
problem.sim.inputU              = 1:length(problem.inputs.ul);
problem.functions_unscaled      = {@L_unscaled, @E_unscaled, @f_unscaled, ...
                                   @g_unscaled, @avrc, @b_unscaled};
problem.data.functions_unscaled = problem.functions_unscaled;
problem.data.ng_eq              = problem.constraints.ng_eq;
problem.constraintErrorTol      = [problem.constraints.gTol_eq, problem.constraints.gTol_neq, ...
                                   problem.constraints.gTol_eq, problem.constraints.gTol_neq, ...
                                   problem.states.xConstraintTol, problem.states.xConstraintTol, ...
                                   problem.inputs.uConstraintTol, problem.inputs.uConstraintTol];

%% ================= UNSCALED COSTS / CONSTRAINTS =================
function stageCost = L_unscaled(~, ~, u, ~, ~, t, vdat)
    thr = u(:,4);
    % Isotropic 1e-6 * ||u_dir||^2 makes all 3 directional eigenvalues strictly positive during coasting
    u_reg = 1e-6 * (u(:,1).^2 + u(:,2).^2 + u(:,3).^2);
    if ischar(vdat.mode_opt) && strcmp(vdat.mode_opt, 'energy')
        stageCost = 0.5 * (thr.^2) + u_reg + 0*t;
    else
        if ischar(vdat.mode_opt), eps_h = 1e-4; else, eps_h = vdat.mode_opt; end
        stageCost = 10.0 * eps_h * (thr.^2) + u_reg + 0*t;
    end
end

function boundaryCost = E_unscaled(~, xf, ~, ~, ~, ~, tf, vdat)
    if ischar(vdat.mode_opt) && strcmp(vdat.mode_opt, 'energy')
        boundaryCost = 0*tf;
    else
        % Factor 10.0: balance bw state, spherical boundary and multipliers
        boundaryCost = -10.0 * xf(7);
    end
end

function bc = b_unscaled(x0, xf, ~, ~, p, ~, tf, vdat, varargin)
    % Arrival epoch shift dtf = p(1) and Departure epoch shift dt0 = p(1) - (tf - tf_guess)
    dtf  = p(1);
    dt0  = dtf - (tf - vdat.tf_guess_tu);

    % Normalized polynomial arguments strictly centered in [-1, 1]
    tau0 = (dt0 - vdat.dt0_mid) / vdat.dt0_half;
    tauf = (dtf - vdat.dtf_mid) / vdat.dtf_half;

    % 4th-degree Horner polynomial evaluation of Earth departure state at tau0
    PE = vdat.P_earth;
    x_earth_t0 = (((PE(1,:)*tau0 + PE(2,:))*tau0 + PE(3,:))*tau0 + PE(4,:))*tau0 + PE(5,:);

    % 4th-degree Horner polynomial evaluation of Target arrival state at tauf
    PT = vdat.P_target;
    x_target_tf = (((PT(1,:)*tauf + PT(2,:))*tauf + PT(3,:))*tauf + PT(4,:))*tauf + PT(5,:);

    bc = [x0(1:6) - x_earth_t0(:); ...
          xf(1:6) - x_target_tf(:)];

    if ~isempty(varargin) && ~isempty(varargin{1})
        vi = varargin{1};
        if length(vi) == 2
            options = vi{1}; t_segment = vi{2};
            if (strcmp(options.discretization,'hpLGR') || ...
                strcmp(options.discretization,'globalLGR')) && options.adaptseg==1
                if size(t_segment,1) > size(t_segment,2)
                    bc = [bc; diff(t_segment)];
                else
                    bc = [bc, diff(t_segment)'];
                end
            end
        end
    end
end

%% ================= SCALED WRAPPERS =================
function stageCost = L(x, xr, u, ur, p, t, vdat)
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
        if ~isempty(xr); xr = scale_variables_back(xr, vdat.Xscale, vdat.Xshift); end
        if ~isempty(ur); ur = scale_variables_back(ur, vdat.Uscale, vdat.Ushift); end
        if isfield(vdat,'Pscale') && ~isempty(p); p = scale_variables_back(p, vdat.Pscale, vdat.Pshift); end
    end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        stageCost = 0*t;
    else
        stageCost = L_unscaled(x, xr, u, ur, p, t, vdat);
    end
end

function boundaryCost = E(x0, xf, u0, uf, p, t0, tf, vdat)
    if isfield(vdat,'Xscale')
        x0 = scale_variables_back(x0', vdat.Xscale, vdat.Xshift);
        xf = scale_variables_back(xf', vdat.Xscale, vdat.Xshift);
        u0 = scale_variables_back(u0', vdat.Uscale, vdat.Ushift);
        uf = scale_variables_back(uf', vdat.Uscale, vdat.Ushift);
        if isfield(vdat,'Pscale') && ~isempty(p); p = scale_variables_back(p', vdat.Pscale, vdat.Pshift)'; end
    end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        boundaryCost = sum(sum(p(:,end-vdat.mode.np*2+1:end)));
    else
        boundaryCost = E_unscaled(x0, xf, u0, uf, p, t0, tf, vdat);
    end
end

function dx = f(x, u, p, t, vdat)
    f_fn = vdat.InternalDynamics;
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
        if isfield(vdat,'Pscale') && ~isempty(p); p = scale_variables_back(p, vdat.Pscale, vdat.Pshift); end
    end
    dx = f_fn(x, u, p, t, vdat);
    if isfield(vdat,'Xscale')
        dx = scale_variables(dx, vdat.Xscale, 0);
    end
end

function c = g(x, u, p, t, vdat)
    g_fn = vdat.InternalDynamics;
    ng   = nargout(g_fn);
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
        if isfield(vdat,'Pscale') && ~isempty(p); p = scale_variables_back(p, vdat.Pscale, vdat.Pshift); end
    end
    if ng == 1
        c = [];
    elseif ng == 2
        [~, c] = g_fn(x, u, p, t, vdat);
    else
        [~, ceq, cneq] = g_fn(x, u, p, t, vdat); c = [ceq, cneq];
    end
    if isfield(vdat,'gFilter'); c(:,vdat.gFilter) = []; end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        c = [c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np), ...
             c+p(:,end-vdat.mode.np+1:end)];
    end
end

function [dx, c] = fg(x, u, p, t, vdat)
    fg_fn = vdat.InternalDynamics;
    ng    = nargout(fg_fn);
    if isfield(vdat,'Xscale')
        x = scale_variables_back(x, vdat.Xscale, vdat.Xshift);
        u = scale_variables_back(u, vdat.Uscale, vdat.Ushift);
        if isfield(vdat,'Pscale') && ~isempty(p); p = scale_variables_back(p, vdat.Pscale, vdat.Pshift); end
    end
    if ng == 1
        c = []; dx = fg_fn(x, u, p, t, vdat);
    elseif ng == 2
        [dx, c] = fg_fn(x, u, p, t, vdat);
    else
        [dx, ceq, cneq] = fg_fn(x, u, p, t, vdat); c = [ceq, cneq];
    end
    if isfield(vdat,'Xscale')
        dx = scale_variables(dx, vdat.Xscale, 0);
    end
    if isfield(vdat,'gFilter'); c(:,vdat.gFilter) = []; end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        c = [c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np), ...
             c+p(:,end-vdat.mode.np+1:end)];
    end
end

function cr = avrc(x, u, p, t, data)
    [cr] = addRateConstraint(x, u, p, t, data);
end

function bc = b(x0, xf, u0, uf, p, t0, tf, vdat, varargin)
    bc = b_unscaled(x0, xf, u0, uf, p, t0, tf, vdat, varargin);
    if isfield(vdat,'Xscale') && ~isempty(bc)
        x0 = scale_variables_back(x0', vdat.Xscale, vdat.Xshift);
        xf = scale_variables_back(xf', vdat.Xscale, vdat.Xshift);
        u0 = scale_variables_back(u0', vdat.Uscale, vdat.Ushift);
        uf = scale_variables_back(uf', vdat.Uscale, vdat.Ushift);
        if isfield(vdat,'Pscale') && ~isempty(p); p = scale_variables_back(p', vdat.Pscale, vdat.Pshift)'; end
        bc = b_unscaled(x0, xf, u0, uf, p, t0, tf, vdat, varargin);
    end
end

function dx = f_unscaled(x, u, p, t, vdat)
    Dyn = vdat.InternalDynamics;
    dx  = Dyn(x, u, p, t, vdat);
end

function c = g_unscaled(x, u, p, t, vdat)
    Dyn = vdat.InternalDynamics;
    ng  = nargout(Dyn);
    if ng == 1
        c = [];
    elseif ng == 2
        [~, c] = Dyn(x, u, p, t, vdat);
    else
        [~, ceq, cneq] = Dyn(x, u, p, t, vdat); c = [ceq, cneq];
    end
end

end