function [x0_guess_nd, xf_nd, tf_bounds_tu, ephem_out] = get_boundary_conditions( ...
    file_earth, file_target, tf_epoch, TOF_guess_days, TOF_window_days, n_rev, units_flag)
%% ========================================================================
%% PROJECT : Interplanetary 3D Low-Thrust Transfer Optimization (ICLOCS2 + IPOPT)
%  AUTHOR  : Gabriele Borghini - November 2026
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
%  Copyright (c) 2026 Gabriele Borghini. All Rights Reserved.
% =========================================================================

%% --- DATA EXTRACTION ---
if nargin < 6, n_rev = []; end
if nargin < 7 || isempty(units_flag), units_flag = 'km_s'; end

AU     = 1.495978707e8;          % [km]
mu_sun = 1.32712440018e11;       % [km^3/s^2]
VU     = sqrt(mu_sun / AU);      % [km/s]
TU     = AU / VU;                % [s]

if ischar(tf_epoch) || isstring(tf_epoch) || isdatetime(tf_epoch)
    JD_tf = juliandate(datetime(tf_epoch));
else
    JD_tf = double(tf_epoch);
end

TOF_min_days = TOF_window_days(1);
TOF_max_days = TOF_window_days(2);

% Departure Julian Dates corresponding to the TOF window: JD_t0 = JD_tf - TOF
JD_t0_guess = JD_tf - TOF_guess_days;
JD_t0_min   = JD_tf - TOF_max_days; % Longest TOF -> Earliest departure
JD_t0_max   = JD_tf - TOF_min_days; % Shortest TOF -> Latest departure

tf_bounds_tu.guess = (TOF_guess_days * 86400) / TU;
tf_bounds_tu.min   = (TOF_min_days   * 86400) / TU;
tf_bounds_tu.max   = (TOF_max_days   * 86400) / TU;

%% --- EPHEMERIS FILE READING ---
data_E = read_ephem_file(file_earth);
data_T = read_ephem_file(file_target);

if JD_t0_min < min(data_E(:,1)) || JD_t0_max > max(data_E(:,1))
    error('Departure window [JD %.2f, %.2f] is out of the Earth ephemeris file range [%.2f, %.2f].', ...
          JD_t0_min, JD_t0_max, min(data_E(:,1)), max(data_E(:,1)));
end
if JD_tf < min(data_T(:,1)) || JD_tf > max(data_T(:,1))
    error('Arrival epoch tf (JD %.2f) is out of the Target ephemeris file range [%.2f, %.2f].', ...
          JD_tf, min(data_T(:,1)), max(data_T(:,1)));
end

%% --- TARGET FINAL STATE AT FIXED ARRIVAL EPOCH (JD_tf) ---
rv_Tf = interp1(data_T(:,1), data_T(:,2:7), JD_tf, 'spline');
if strcmpi(units_flag, 'au_d')
    r_Tf_km = rv_Tf(1:3) * AU;  v_Tf_kms = rv_Tf(4:6) * (AU / 86400);
else
    r_Tf_km = rv_Tf(1:3);       v_Tf_kms = rv_Tf(4:6);
end
sphf_raw = cart2sph_dynamics(r_Tf_km / AU, v_Tf_kms / VU);

%% --- SMOOTH C^inf ANALYTICAL MODEL OF EARTH DEPARTURE AS FUNCTION OF TOF ---
N_sample  = 200;
t_tu_grid = linspace(tf_bounds_tu.min, tf_bounds_tu.max, N_sample).';
% For a given duration t_tu, departure epoch is JD_tf - t_days
JD_dep_grid = JD_tf - (t_tu_grid * TU) / 86400;

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

sph0_guess  = interp1(t_tu_grid, sph_E_grid, tf_bounds_tu.guess, 'spline');
dtheta_base = mod(sphf_raw(2) - sph0_guess(2), 2*pi);

r_mean     = 0.5 * (sph0_guess(1) + sphf_raw(1));
dtheta_est = (1.0 / (r_mean^1.5)) * tf_bounds_tu.guess;
n_rev_best = max(0, round((dtheta_est - dtheta_base) / (2*pi)));

if isempty(n_rev)
    n_rev = n_rev_best;
elseif n_rev ~= n_rev_best
    warning('n_rev=%d enforced, automatic estimate would give %d.', n_rev, n_rev_best);
end

% Align final target longitude consistently with Earth departure guess + n_rev
sphf    = sphf_raw;
sphf(2) = sph0_guess(2) + dtheta_base + n_rev * 2*pi;

% 4th-degree polynomial fit of Earth initial state centered at tf_guess
dt_grid     = t_tu_grid - tf_bounds_tu.guess;
poly_coeffs = zeros(5, 6);
for k = 1:6
    poly_coeffs(:,k) = polyfit(dt_grid, sph_E_grid(:,k), 4).';
end

pp_earth    = spline(t_tu_grid, sph_E_grid.');
sph0_fit    = poly_coeffs(5,:);
x0_guess_nd = [sph0_fit, 1.0];
xf_nd       = [sphf, 1.0];

%% --- DIAGNOSTIC OUTPUT ---
ephem_out.JD_tf         = JD_tf;
ephem_out.JD_t0_guess   = JD_t0_guess;
ephem_out.pp_earth      = pp_earth;
ephem_out.poly_coeffs   = poly_coeffs;
ephem_out.tf_guess_tu   = tf_bounds_tu.guess;
ephem_out.dtheta_deg    = rad2deg(xf_nd(2) - x0_guess_nd(2));
ephem_out.TU            = TU;
ephem_out.n_rev         = n_rev;
ephem_out.n_rev_est     = n_rev_best;
ephem_out.dtheta_margin = (dtheta_est - dtheta_base) / (2*pi) - n_rev;
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