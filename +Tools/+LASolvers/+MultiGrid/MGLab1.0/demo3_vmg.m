%==============================================
% Multigrid V-Cycle Solver, adapted for Demo 3.
%==============================================

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x,resids,its,pde_error,rho] = demo3_vmg(A,b,x,xtrue,rtol,prtol,max_it,max_time,max_mflop)

r0 = b - A * x;    r = r0;
rn0 = norm(r0);    rn = rn0;

iter = 0;
results=update_results([],'VMG',iter,rn);

pde_error = [norm(xtrue - x)];

while (rn > rtol * rn0 & iter < max_it)
    level = 1;
    x = vmg_cycle(level,b,x);
    r = b - A*x;
    rn = norm(r); 
    iter = iter + 1;
    results=update_results(results,'VMG',iter,rn);
    pde_error = [pde_error; norm(xtrue-x)];
end

its=results(:,1);      resids=results(:,4);

rho=10^(log10(resids(10)/resids(2))/9);

fprintf('rho = %g.\n', rho)
