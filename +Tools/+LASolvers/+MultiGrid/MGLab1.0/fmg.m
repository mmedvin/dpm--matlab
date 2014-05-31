%FMG    Full Multigrid Solver
%
%       [X,RESIDS,ITS] = FMG(A,B) uses the full-multigrid cycle to 
%       solve the linear system AX=B.  RESIDS is the final residual
%       and ITS is 1.
%
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x,resids,its] = fmg(A,b,x)

x = fmg_cycle(1,b);
rn = norm(b-A*x);
results=update_results([],'Full Multigrid',1,rn);

its=results(:,1);
resids=results(:,4);

fprintf('Residual reduction = %g.\n', rn)
