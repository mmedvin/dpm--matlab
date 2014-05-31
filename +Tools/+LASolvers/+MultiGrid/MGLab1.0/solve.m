%SOLVE  Solve a linear system.
%
%       [X,RESIDS,ITS] = SOLVE(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP,...
%       RESTART) applies a solver defined by "solver_flag", with the given
%       tolerances and limits, to a linear system AX=B.  The solution X, 
%       residual history RESIDS, and iterations ITS are returned.
%
%       Accesses global variables in "include_flags"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x,resids,its] = solve(A,b,x0,...
    rtol,prtol,max_it,max_time,max_mflop,restart)

include_flags 

disp (sprintf('Running...\n'));
if solver_flag == VMG
    [x,resids,its] = vmg (A,b,x0,rtol,prtol,max_it,max_time,max_mflop);
elseif solver_flag == FMG
    [x,resids,its] = fmg (A,b);
elseif solver_flag == PCG
    [x,resids,its] = pcg (A,b,x0,rtol,prtol,max_it,max_time,max_mflop);
elseif solver_flag == BICG_STAB
    [x,resids,its] = pbicgstab (A,b,x0,rtol,prtol,max_it,max_time,max_mflop);
elseif solver_flag == CGS
    [x,resids,its] = pcgs (A,b,x0,rtol,prtol,max_it,max_time,max_mflop);
elseif solver_flag == GMRES
    [x,resids,its] = pgmres (A,b,x0,rtol,prtol,max_it,max_time,max_mflop,restart);
elseif solver_flag == SOR
    [x,resids,its] = sor (A,b,x0,rtol,prtol,max_it,max_time,max_mflop);
end

fprintf('Relative residual = %g \n', norm(b-A*x))
disp (sprintf('Done.\n'));
