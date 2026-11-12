function [x0_guess_nd, xf_guess_nd, tf_bounds_tu, dtf_bounds_tu, ephem_out] = get_boundary_conditions( ...
    file_earth, file_target, tf_epoch, arr_window_days, TOF_guess_days, TOF_window_days, n_rev, units_flag)
%% ========================================================================
%% PROJECT : Interplanetary 3D Low-Thrust Transfer Optimization (ICLOCS2 + IPOPT)
%  AUTHOR  : Gabriele Borghini - November 2026
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
%  Ephemeris parser and boundary condition generator for fixed-arrival / 
%  free-departure (free-TOF) transfers.
%   - Parses JPL Horizons vector tables (.txt/.csv) for Earth and target bodies
%   - Evaluates the fixed target state at the prescribed arrival epoch (tf_epoch)
%   - Automatically aligns true longitude theta and revolution count (n_rev)
%   - Generates smooth C^inf local polynomial approximations for Earth departure matching
%
%  Inputs:   file_earth, file_target, tf_epoch, TOF_guess_days, TOF_window_days, n_rev
%  Outputs:  x0_guess_nd, xf_nd, tf_bounds_tu, ephem_out
% =========================================================================

%% --- DATA EXTRACTION ---
if nargin < 7, n_rev = []; end
if nargin < 8 || isempty(units_flag), units_flag = 'km_s'; end

AU     = 1.495978707e8;          % [km]
mu_sun = 1.32712440018e11;       % [km^3/s^2]
VU     = sqrt(mu_sun / AU);      % [km/s]
TU     = AU / VU;                % [s]

if ischar(tf_epoch) || isstring(tf_epoch) || isdatetime(tf_epoch)
    JD_tf_guess = juliandate(datetime(tf_epoch));
else
    JD_tf_guess = double(tf_epoch);
end

dtf_min_days = arr_window_days(1);
dtf_max_days = arr_window_days(2);

TOF_min_days = TOF_window_days(1);
TOF_max_days = TOF_window_days(2);

% Arrival Julian Dates corresponding to the arrival window around JD_tf_guess
JD_tf_min = JD_tf_guess + dtf_min_days;
JD_tf_max = JD_tf_guess + dtf_max_days;

% Departure Julian Dates: JD_t0 = JD_tf - TOF
JD_t0_guess = JD_tf_guess - TOF_guess_days;
JD_t0_min   = JD_tf_min   - TOF_max_days;
JD_t0_max   = JD_tf_max   - TOF_min_days;

% Bounds on Time of Flight tf [TU]
tf_bounds_tu.guess = (TOF_guess_days * 86400) / TU;
tf_bounds_tu.min   = (TOF_min_days   * 86400) / TU;
tf_bounds_tu.max   = (TOF_max_days   * 86400) / TU;

% Bounds on Arrival Epoch Shift p(1) = dtf [TU] relative to JD_tf_guess
dtf_bounds_tu.min = (dtf_min_days * 86400) / TU;
dtf_bounds_tu.max = (dtf_max_days * 86400) / TU;
% Keep guess near 0.0 (nominal arrival date), just 10% inside the bound if 0.0 is an endpoint
margin_tu = 0.1 * (dtf_bounds_tu.max - dtf_bounds_tu.min);
dtf_bounds_tu.guess = min(max(0.0, dtf_bounds_tu.min + margin_tu), dtf_bounds_tu.max - margin_tu);

% Resulting Departure Shift range dt0 [TU] relative to JD_t0_guess:
dt0_min_tu = dtf_bounds_tu.min - (tf_bounds_tu.max - tf_bounds_tu.guess);
dt0_max_tu = dtf_bounds_tu.max - (tf_bounds_tu.min - tf_bounds_tu.guess);

%% --- EPHEMERIS FILE READING ---
data_E = read_ephem_file(file_earth);
data_T = read_ephem_file(file_target);

if JD_t0_min < min(data_E(:,1)) || JD_t0_max > max(data_E(:,1))
    error('Departure window [JD %.2f, %.2f] is out of the Earth ephemeris file range [%.2f, %.2f].', ...
          JD_t0_min, JD_t0_max, min(data_E(:,1)), max(data_E(:,1)));
end
if JD_tf_min < min(data_T(:,1)) || JD_tf_max > max(data_T(:,1))
    error('Arrival window [JD %.2f, %.2f] is out of the Target ephemeris file range [%.2f, %.2f].', ...
          JD_tf_min, JD_tf_max, min(data_T(:,1)), max(data_T(:,1)));
end

%% --- SMOOTH C^inf ANALYTICAL MODEL OF EARTH DEPARTURE AS FUNCTION OF dt0 ---
N_sample    = 200;
dt0_tu_grid = linspace(dt0_min_tu, dt0_max_tu, N_sample).';
JD_dep_grid = JD_t0_guess + (dt0_tu_grid * TU) / 86400;

rv_E_grid = interp1(data_E(:,1), data_E(:,2:7), JD_dep_grid, 'spline');
if strcmpi(units_flag, 'au_d')
    r_E_nd = rv_E_grid(:,1:3);
    v_E_nd = rv_E_grid(:,4:6) * (AU / 86400) / VU;
else
    r_E_nd = rv_E_grid(:,1:3) / AU;
    v_E_nd = rv_E_grid(:,4:6) / VU;
end

sph_E_grid = zeros(N_sample, 6);
for i = 1:N_sample
    sph_E_grid(i,:) = cart2sph_dynamics(r_E_nd(i,:), v_E_nd(i,:));
end
sph_E_grid(:,2) = unwrap(sph_E_grid(:,2));
sph_E_grid(:,6) = unwrap(sph_E_grid(:,6));

sph0_guess = interp1(dt0_tu_grid, sph_E_grid, dtf_bounds_tu.guess, 'spline');

%% --- SMOOTH C^inf ANALYTICAL MODEL OF TARGET ARRIVAL AS FUNCTION OF dtf ---
dtf_tu_grid = linspace(dtf_bounds_tu.min, dtf_bounds_tu.max, N_sample).';
JD_arr_grid = JD_tf_guess + (dtf_tu_grid * TU) / 86400;

rv_T_grid = interp1(data_T(:,1), data_T(:,2:7), JD_arr_grid, 'spline');
if strcmpi(units_flag, 'au_d')
    r_T_nd = rv_T_grid(:,1:3);
    v_T_nd = rv_T_grid(:,4:6) * (AU / 86400) / VU;
else
    r_T_nd = rv_T_grid(:,1:3) / AU;
    v_T_nd = rv_T_grid(:,4:6) / VU;
end

sph_T_grid = zeros(N_sample, 6);
for i = 1:N_sample
    sph_T_grid(i,:) = cart2sph_dynamics(r_T_nd(i,:), v_T_nd(i,:));
end
sph_T_grid(:,2) = unwrap(sph_T_grid(:,2));
sph_T_grid(:,6) = unwrap(sph_T_grid(:,6));

sphf_raw    = interp1(dtf_tu_grid, sph_T_grid, dtf_bounds_tu.guess, 'spline');
dtheta_base = mod(sphf_raw(2) - sph0_guess(2), 2*pi);

r_mean     = 0.5 * (sph0_guess(1) + sphf_raw(1));
dtheta_est = (1.0 / (r_mean^1.5)) * tf_bounds_tu.guess;
n_rev_best = max(0, round((dtheta_est - dtheta_base) / (2*pi)));

if isempty(n_rev)
    n_rev = n_rev_best;
elseif n_rev ~= n_rev_best
    warning('n_rev=%d enforced, automatic estimate would give %d.', n_rev, n_rev_best);
end

% Align target longitude curve consistently with Earth departure guess + n_rev
lon_offset      = (sph0_guess(2) + dtheta_base + n_rev * 2*pi) - sphf_raw(2);
sph_T_grid(:,2) = sph_T_grid(:,2) + lon_offset;

% Centered and scaled 4th-degree polynomial fits strictly in [-1, 1] (fast Hessian + high accuracy)
dt0_mid  = 0.5 * (dt0_min_tu + dt0_max_tu);
dt0_half = max(0.5 * (dt0_max_tu - dt0_min_tu), 1e-6);

dtf_mid  = 0.5 * (dtf_bounds_tu.min + dtf_bounds_tu.max);
dtf_half = max(0.5 * (dtf_bounds_tu.max - dtf_bounds_tu.min), 1e-6);

tau_E_grid = (dt0_tu_grid - dt0_mid) / dt0_half;
tau_T_grid = (dtf_tu_grid - dtf_mid) / dtf_half;

poly_earth  = zeros(5, 6);
poly_target = zeros(5, 6);
for k = 1:6
    poly_earth(:,k)  = polyfit(tau_E_grid, sph_E_grid(:,k), 4).';
    poly_target(:,k) = polyfit(tau_T_grid, sph_T_grid(:,k), 4).';
end

pp_earth  = spline(dt0_tu_grid, sph_E_grid.');
pp_target = spline(dtf_tu_grid, sph_T_grid.');

x0_guess_nd = [ppval(pp_earth,  dtf_bounds_tu.guess).', 1.0];
xf_guess_nd = [ppval(pp_target, dtf_bounds_tu.guess).', 1.0];

%% --- DIAGNOSTIC OUTPUT ---
ephem_out.JD_tf_guess     = JD_tf_guess;
ephem_out.JD_t0_guess     = JD_t0_guess;
ephem_out.JD_tf_window    = [JD_tf_min, JD_tf_max];
ephem_out.JD_t0_window    = [JD_t0_min, JD_t0_max];
ephem_out.dep_window_days = [JD_t0_min - JD_t0_guess, JD_t0_max - JD_t0_guess];
ephem_out.pp_earth        = pp_earth;
ephem_out.pp_target       = pp_target;
ephem_out.poly_earth      = poly_earth;
ephem_out.poly_target     = poly_target;
ephem_out.dt0_mid         = dt0_mid;
ephem_out.dt0_half        = dt0_half;
ephem_out.dtf_mid         = dtf_mid;
ephem_out.dtf_half        = dtf_half;
ephem_out.tf_guess_tu     = tf_bounds_tu.guess;
ephem_out.dtheta_deg      = rad2deg(xf_guess_nd(2) - x0_guess_nd(2));
ephem_out.TU              = TU;
ephem_out.n_rev           = n_rev;
ephem_out.n_rev_est       = n_rev_best;
ephem_out.dtheta_margin   = (dtheta_est - dtheta_base) / (2*pi) - n_rev;
end

function sph = cart2sph_dynamics(r_vec, v_vec)
    r   = norm(r_vec);
    v   = norm(v_vec);
    phi   = asin(max(-1.0, min(1.0, r_vec(3) / r)));
    theta = atan2(r_vec(2), r_vec(1));
    
    u_U = r_vec(:) / r;
    u_E = [-sin(theta); cos(theta); 0];
    u_N = [-sin(phi)*cos(theta); -sin(phi)*sin(theta); cos(phi)];
    
    v_U = dot(v_vec(:), u_U);
    v_N = dot(v_vec(:), u_N);
    v_E = dot(v_vec(:), u_E);
    
    gamma = asin(max(-1.0, min(1.0, v_U / v)));
    chi   = atan2(v_E, v_N);
    sph   = [r, theta, phi, v, gamma, chi];
end

function data = read_ephem_file(filepath)
    txt = fileread(filepath);
    soe_idx = strfind(txt, '$$SOE');
    eoe_idx = strfind(txt, '$$EOE');
    if ~isempty(soe_idx) && ~isempty(eoe_idx)
        block = txt(soe_idx(1)+5 : eoe_idx(1)-1);
        lines = splitlines(strtrim(block));
        data  = nan(length(lines), 7);
        for i = 1:length(lines)
            line_str = strtrim(lines{i});
            if isempty(line_str), continue; end
            parts = strsplit(line_str, ',');
            if length(parts) >= 8
                data(i,:) = [str2double(parts{1}), str2double(parts(3:8))];
            else
                tokens = regexp(line_str, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
                nums   = str2double(tokens);
                if length(nums) >= 7
                    data(i,:) = [nums(1), nums(end-5:end)];
                end
            end
        end
        data = data(~any(isnan(data), 2), :);
    else
        raw  = readmatrix(filepath);
        raw  = raw(~any(isnan(raw), 2), :);
        data = [raw(:,1), raw(:, end-5:end)];
    end
end