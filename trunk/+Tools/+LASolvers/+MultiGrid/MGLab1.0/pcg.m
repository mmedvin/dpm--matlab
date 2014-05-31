%PCG    Preconditioned conjugate gradient method
%
%       [X,RESIDS,ITS]=PCG(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP)
%       solves the system AX = B using the preconditioned conjugate gradient 
%       method with the given tolerances and limits.  A should be symmetric
%       positive definite.
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x,resids,its] = pcg(A,b,x0,rtol,prtol,max_it,max_time,max_mflop)

x = x0; 
r = b - A*x; 
z = precondition(A,r);
prn = norm(z);

iter = 0;
results=update_results([],'PCG',iter,prn);

bn = 1;
rn = 1;
pbn = precondition(A,b);

% z could be named pr but we use the more common notation.

while (~converged(bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop))

   rho = r' * z;
   if iter == 0
      p = z;
   else
      beta = rho / rho_old;
      p = z + beta * p;
   end
   rho_old = rho;

   Ap = A*p;
   mu = p' * Ap;
   alpha = rho / mu;
   x = x + alpha *  p;
   r = r - alpha * Ap;
   z = precondition (A,r);

   prn = norm(z);

   iter = iter+1;
   results=update_results(results,'PCG',iter,prn);
end

its=results(:,1);
resids=results(:,4);
