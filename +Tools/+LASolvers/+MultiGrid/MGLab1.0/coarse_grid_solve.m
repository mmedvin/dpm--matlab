%COARSE_GRID_SOLVE Solve the coarse-grid system
%
%       COARSE_GRID_SOLVE(LEVEL,B) solves the linear system at the grid
%       level LEVEL with the right-hand side B.  How the system is solved
%       depends on the global variable "coarse_solver_flag" as follows:
%
%          DIRECT        - Sparse Gaussian elimination.
%          SMOOTHER      - A constant number of applications of the smoother
%          PCG           - The PCG method (NOT IMPLEMENTED).
%          BICG_STAB     - The PBICG_STAB method (NOT IMPLEMENTED).
%          GMRES         - The PGMRES method (NOT IMPLEMENTED).
%          
%       Accesses global variables in "include_flags"
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function u_out = coarse_grid_solve(level, b)

include_flags
include_globals

if (coarse_solver_flag == DIRECT)
    eval(['u_out = A', num2str(level), ' \ b;']);
elseif (coarse_solver_flag == SMOOTHER)
    u_out  = smooth(level, b, b, 'coarse');
elseif (coarse_solver_flag == PCG)
    disp(sprintf('PCG coarse-grid solve not implemented, using DIRECT'));
    eval(['u_out = A', num2str(level), ' \ b;']);
elseif (coarse_solver_flag == BICG_STAB)
    disp(sprintf('BICG_STAB coarse-grid solve not implemented, using DIRECT'));
    eval(['u_out = A', num2str(level), ' \ b;']);
elseif (coarse_solver_flag == GMRES)
    disp(sprintf('GMRES coarse-grid solve not implemented, using DIRECT'));
    eval(['u_out = A', num2str(level), ' \ b;']);
end
