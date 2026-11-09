function options = settings_low_thrust(varargin)
%% ========================================================================
%% PROJECT : Interplanetary 3D Low-Thrust Transfer Optimization (ICLOCS2 + IPOPT)
%  AUTHOR  : Gabriele Borghini - November 2026
%  Based on ICLOCS2 settings template (Copyright (c) 2018 Imperial College London, MIT License)
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
%  Configuration parameters for ICLOCS2 direct collocation (trapezoidal transcription) 
%  and IPOPT optimization engine tailored for low-thrust heliocentric transfers.
%
%  Usage:
%   options = problem.settings(300)    -> trapezoidal method, 300 nodes
% =========================================================================
%  Copyright (c) 2026 Gabriele Borghini. All Rights Reserved.
% =========================================================================

%% Transcription Method
options.transcription    = 'direct_collocation';
options.min_res_mode     = 'alternating';
options.min_res_priority = 'low_res_error';
options.resCostWeight    = [];
options.errortype        = 'local_abs';

%% Discretization Method
% 'trapezoidal' eliminates hidden internal variables at Hermite midpoints,
% making the control profile in solution.U identical to the one integrated by ode45!
options.discretization = 'trapezoidal';
options.resultRep      = 'default';

%% Derivative Generation
options.derivatives      = 'numeric';
% options.derivatives    = 'adigator';
options.adigatorPath     = 'D:\OptimalControl\adigator-master';
options.perturbation.H   = [];
options.perturbation.J   = [];

%% NLP Solver
options.NLPsolver = 'ipopt';

% IPOPT settings
options.ipopt.tol                         = 1e-7;
options.ipopt.constr_viol_tol             = 1e-7;
options.ipopt.acceptable_tol              = 1e-6;
options.ipopt.acceptable_constr_viol_tol  = 1e-6;
options.ipopt.acceptable_iter             = 15;
options.ipopt.max_iter                    = 1000;
options.ipopt.print_level                 = 5;
options.ipopt.mu_strategy                 = 'adaptive';
options.ipopt.hessian_approximation       = 'exact';
options.ipopt.limited_memory_max_history  = 6;
options.ipopt.limited_memory_max_skipping = 1;

% fmincon settings (fallback)
options.fmincon.optimoptions = optimoptions('fmincon', 'Display', 'iter', ...
    'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
    'MaxIterations', 5000);
options.fmincon.hessian_approximation = 'exact';

%% Meshing Strategy
options.meshstrategy    = 'fixed';
options.MeshRefinement  = 'Auto';
options.MRstrategy      = 'aggressive';
options.maxMRiter       = 3;
options.disContInputs   = 0;
options.mintimeinterval = 1e-4;
options.maxtimeinterval = inf;
options.tau             = 0;

%% Other Settings
options.start           = 'Cold';
options.scaling         = 0;
options.reorderLGR      = 0;
options.resminEarlyStop = 0;
options.ECH.enabled     = 0;
options.ECH.buffer_pct  = 0.1;
options.regstrategy     = 'off';

%% Auto h/hp selection (DO NOT MODIFY - from original ICLOCS2 template)
if nargin == 2
    if strcmp(varargin{2}, 'h')
        options.nodes          = varargin{1};
        options.discretization = 'trapezoidal';
    else
        if length(varargin{1}) == 1
            options.nsegment = varargin{1};
            options.pdegree  = varargin{2};
        else
            options.tau_segment = varargin{1};
            options.npsegment   = varargin{2};
        end
        options.discretization = 'hpLGR';
    end
else
    options.nodes = varargin{1};
end

%% Output Settings
options.print.time                 = 1;
options.print.relative_local_error = 1;
options.print.cost                 = 1;
options.plot                       = 0;
end