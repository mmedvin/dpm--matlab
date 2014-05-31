%RUN    Apply a solver to a linear system.
%
%       [X1,RESIDS1,ITS1] = RUN applies a solver to a discretized PDE problem
%       and returns the solution X1, history of residuals RESIDS1, and
%       iteration vector ITS1.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x1,resids1,its1] = run

include_globals
include_flags

show_params; drawnow;

if (generate_matrix)
   [A1,N1] = get_matrix(nx1,ny1);
   generate_matrix = 0;
   multigrid_setup;
end

if (generate_rhs)
   b1 = get_rhs(nx1,ny1);   
   generate_rhs = 0;
end

x0 = zeros(N1,1);
[x1,resids1,its1] = solve (A1,b1,x0,rtol,prtol,max_it,max_time,max_mflop,restart);
