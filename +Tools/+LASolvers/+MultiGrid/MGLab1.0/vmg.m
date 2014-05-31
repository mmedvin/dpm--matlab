%VMG    Multigrid V-Cycle Solver
%
%       [X,RESIDS,ITS]=VMG(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP)
%       solves the system AX = B iteratively using multigrid cycles whose type
%       is defined by "cycle_flag".  The stopping criteria is given by
%       the input tolerances and limits.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x,resids,its] = vmg(A,b,x,rtol,prtol,max_it,max_time,max_mflop)

r0 = b - A * x; 
r = r0;
rn0 = norm(r0);  
rn = rn0;

iter = 0;
results=update_results([],'V-Cycle',iter,rn);

bn=norm(b);
pbn = bn;
prn = rn;

% We stop on the true resid, but the test is actually performed on the precond
% resid.

while (~converged(bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop))
    level = 1;
    x = mg_cycle(level,b,x);
    r = b - A*x;
    iter = iter + 1;
    rn = norm(r);
    prn = rn;
    results=update_results(results,'V-Cycle',iter,prn);
end

its=results(:,1);
resids=results(:,4);

rho=10^(log10(resids(iter+1)/resids(2))/((iter+1)-2));
fprintf('rho = %g.\n', rho)
