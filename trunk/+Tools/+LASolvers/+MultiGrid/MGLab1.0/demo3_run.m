%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function demo3_run

demo_globals
include_globals
include_flags

subplot(1,1,1), cla, hold off

if DEMO_VAR1 == 1
   NX1 = [7 15];
elseif DEMO_VAR1 == 2
   NX1 = [7 15 31];
elseif DEMO_VAR1 == 3
   NX1 = [7 15 31 63];
else
   disp('demo3_run')
end
for nx1 = NX1
    ny1 =nx1;                    % Set fine mesh size


    % Problem parameter defaults

    prob_args = [0 0];

    % Mesh parameter defaults

    coarse_level = 3;            % Set current maximum number of levels

    % Smoother parameter defaults

    nu1 = 1; nu2 = 1;            % Set number of pre- and post-smoothings
    wt = 0.95;                   % Set smoother weighting

    % Selection flag defaults

    problem_flag = POISSON;      % Select continuous problem
    rhs_flag = 1;                % Select right-hand side
    smooth_flag = GAUSS_SEIDEL;  % Select smoother
    restrict_flag = HALF_WEIGHTING; % Select restriction operator
    interp_flag = CUBIC;         % Select interpolation operator
    solver_flag = VMG;           % Select solver
    precon_flag = NONE;          % Select preconditioner
    cycle_flag = V_CYCLE;        % Select cycle flag
    coarsening_flag = 0;         % Select coarsening
    coarse_solver_flag = DIRECT; % Select coarse-grid solve 
    cycle_flag = V_CYCLE;        % Select MG cycle

    % Method parameter defaults

    rtol = 1e-16;           % Set stopping tolerance on weighted residual
    prtol = 1e-16;          % Set stopping tol. on weighted pseudo-resid.
    max_it = 12;                 % Set maximum number of iterations
    max_time = 0;                % Set maximum number of seconds
    max_mflop = 0;               % Set maximum number of mflops
    num_runs = 1;                % Set number of experiments to run

    % Linear System defaults

    generate_matrix = 1;         % Must generate the matrix initially
    generate_rhs = 1;        % Must generate the right-hand side initially
    matrix_type = [];            % Matrix properties initially unknown


    if (generate_matrix)
       [A1,N1] = get_matrix(nx1,ny1);
       generate_matrix = 0;
       multigrid_setup;
    end
    
    if (generate_rhs)
       h = 1 / (nx1+1);
       [XX,YY] = meshgrid([h:h:(1-h)]);
       UTRUE = sin(2*pi*XX) .* sin(3*pi*YY);
       utrue = reshape(UTRUE,N1,1);
       RHS = (h*h) * 13*pi*pi*UTRUE;
       b = reshape(RHS,N1,1);
       generate_rhs = 0;
    end

    x0 = zeros(N1,1);
    [x,resids,its,pde_error,rho] = demo3_vmg(A1,b,x0,utrue,...
			  rtol,prtol,max_it,max_time,max_mflop);
    
    xx = [0:length(resids)-1]';

    if nx1 == NX1(1)
       subplot(221)
    elseif nx1 == NX1(2)
       subplot(222)
    elseif nx1 == NX1(3)
       subplot(223)
    elseif nx1 == NX1(4)
       subplot(224)
    end

    hold off
    semilogy([0,20],[1e1, 1e-18],'.')
    hold on
    semilogy(xx,resids/resids(1),'r')
    semilogy(xx,pde_error/pde_error(1),'b')

    title(['Poisson''s equation on a ', num2str(nx1),' x ',num2str(nx1), ' mesh'])
    text(5,3,['PDE error,  (', num2str(-log10(pde_error(max_it+1)/pde_error(1))), ' digits)'])
    text(10,1e-8,'Discrete residual')
    %text(2, 1e-14, ['Residual reduction factor = ', num2str(rho)])

    %text(10,10, [num2str(-log10(pde_error(max_it+1)/pde_error(1))), ' digits'])
    text(2, 1e-14, num2str(rho))
    hold off
    pause(2)
end


% Reset axes
% close(fig1); close; MGLab
