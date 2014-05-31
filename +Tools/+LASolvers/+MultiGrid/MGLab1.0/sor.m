%SOR    Successive Over-Relaxation method
%
%       [X,RESIDS,ITS]=SOR(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP)
%       solves the system AX = B using the successive over-relaxation
%       method with the given tolerances and limits.  Calls get_SOR_omega
%       which reads omega from globals.
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 20 June 1995

function [x,resids,its] = sor(A,b,x0,rtol,prtol,max_it,max_time,max_mflop)

omega = get_SOR_omega;

[N,N2]=size(A);

D = spdiags(spdiags(A,[0]),[0],N,N);
L = tril(A,-1);
M = (1/omega) * D + L;
clear D, L;

x = x0;
r = b - A * x; 
rn = norm(r);  

iter = 0;
results=update_results([],'SOR',iter,rn);

bn=norm(b);
pbn = bn;
prn = rn;

% We stop on the true resid, but the test is actually performed 
% on the precond residual.

while (~converged(bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop))
    update = M \ r;
    x = x + update; 
    r = r - A * update;
    iter = iter + 1;
    rn = norm(r);
    prn = rn;
    results=update_results(results,'SOR',iter,prn);
end

its=results(:,1);
resids=results(:,4);

rho=10^(log10(resids(iter+1)/resids(2))/((iter+1)-2));
fprintf('rho = %g.\n', rho)
