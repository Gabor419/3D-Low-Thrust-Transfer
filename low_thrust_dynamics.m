function [dx, ceq] = low_thrust_dynamics(x, u, ~, ~, vdat)
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
%  3D Heliocentric Inertial Dynamics (spherical coordinates) for low-thrust 
%  interplanetary transfers.
%  Optimized version to help convergence:
%   - Trigonometric terms evaluated once to streamline the computation graph
%   - Inertial frame formulation (no planetary rotation, atmospheric drag, or staging)
%   - Singularity-free local Cartesian steering controls and variable engine throttling
%
%  State:   x = [r, theta, phi, v, gamma, chi, m]   [AU, rad, rad, VU, rad, rad, MU]
%  Control: u = [ux, uy, uz, throttle]              [-, -, -, -]
%
%  [ux, uy, uz] = local unit thrust direction vector (Up, North, East)
%  throttle     = thrust modulation factor in [0, 1]
% =========================================================================

%% --- CONSTANTS ---
mu       = vdat.mu;        % = 1.0 [AU^3/TU^2]
Tmax     = vdat.Tmax;      % Dimensionless max thrust
mdot_max = vdat.mdot_max;  % Dimensionless max mass flow rate

r     = x(:,1);  phi   = x(:,3);  v   = x(:,4);
gamma = x(:,5);  chi   = x(:,6);  m   = x(:,7);

ux       = u(:,1);
uy       = u(:,2);
uz       = u(:,3);
throttle = u(:,4);

%% --- TRIGONOMETRIC FUNCTIONS ---
sg = sin(gamma);  cg = cos(gamma);
sc = sin(chi);    cc = cos(chi);
sp = sin(phi);    cp = cos(phi);
uh = uy.*cc + uz.*sc;                        % Horizontal projection along velocity azimuth
uc = uz.*cc - uy.*sc;                        % Cross-track horizontal projection

%% --- THRUST AND MASS FLOW RATE ---
T    = Tmax .* throttle;
mdot = mdot_max .* throttle;

T_v     = T .* (ux .* sg + uh .* cg);
T_gamma = T .* (ux .* cg - uh .* sg);
T_chi   = T .* uc;

%% --- EQUATIONS OF MOTION (INERTIAL SPHERICAL REF. FRAME) ---
v2  = v.^2 + 1e-12;
v_s = sqrt(v2);

g   = mu ./ (r.^2);
vcg = v_s .* cg;
vr  = v2 ./ r;

dr     = v_s .* sg;
dtheta = (vcg .* sc) ./ (r .* cp);
dphi   = (vcg .* cc) ./ r;
dv     = (T_v ./ m) - g .* sg;
dgamma = (T_gamma ./ m + (vr - g) .* cg) ./ v_s;
dchi   = (T_chi ./ m + ((vr .* cg) .* cg) .* sc .* (sp ./ cp)) ./ vcg;
dm     = -mdot;

dx  = [dr, dtheta, dphi, dv, dgamma, dchi, dm];
ceq = ux.^2 + uy.^2 + uz.^2 - 1.0;
end