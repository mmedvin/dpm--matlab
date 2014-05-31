%SET_DEFAULTS Initialize global variables
%
%       This M-file initializes global variables in "include_globals",
%       "include_flags", and "include_figs" to their default values.
%
%       Accesses global variables in "include_globals"
%       Accesses global variables in "include_flags"
%       Accesses global variables in "include_figs"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

include_globals
include_flags
include_figs

% Mesh parameter defaults

    nx1 = 15; ny1 =15;           % Set fine mesh size
    coarse_level = 3;            % Set current maximum number of levels

% Problem parameter defaults

    prob_args = [0 0];

% Smoother parameter defaults

    nu1 = 1; nu2 = 1;            % Set number of pre- and post-smoothings
    wt = 0.95;                   % Set smoother weighting

% Selection flag defaults

    problem_flag = POISSON;      % Select continuous problem
    rhs_flag = 1;                % Select right-hand side
    smooth_flag = WEIGHTED_JACOBI;  % Select smoother
    restrict_flag = HALF_WEIGHTING; % Select restriction operator
    interp_flag = CUBIC;         % Select interpolation operator
    solver_flag = VMG;           % Select solver
    precon_flag = MG_CYCLE;       % Select preconditioner
    cycle_flag = V_CYCLE;        % Select MG cycle
    coarsening_flag = 0;         % Select coarsening
    coarse_solver_flag = DIRECT; % Select coarse-grid solve 

% Method parameter defaults

    rtol = 0;                    % Set stopping tolerance on weighted residual
    prtol = 1e-4;                % Set stopping tol. on weighted pseudo-resid.
    max_it = 5;                  % Set maximum number of iterations
    max_time = 0;                % Set maximum number of seconds
    max_mflop = 0;               % Set maximum number of mflops
    num_runs = 1;                % Set number of experiments to run
    restart = 5;                 % Set GMRES restart parameter

% Linear System defaults

    generate_matrix = 1;         % Must generate the matrix initially
    generate_rhs = 1;            % Must generate the right-hand side initially
    matrix_type = [];            % Matrix properties initially unknown

% Figures

    param_fig=0;

    %main_position = [300,30,400,300];
    %main_position = [20,30, 1110, 820];
    %main_position = [ 460   450   570   370];
    %param_position = [ 90   450   360   400];

    main_position = [390 310 750 550];
    param_position = [10 455 370 400];

