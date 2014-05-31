%PBICGSTAB Preconditioned stabilized bi-conjugate gradient method.
%
%       [X,RESIDS,ITS]=PBICGSTAB(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP)
%       solves the system AX = B using the preconditioned stabilized bi-
%       conjugate gradient method with the given tolerances and limits.
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x,resids,its]=pbicgstab(A,b,x0,rtol,prtol,max_it,max_time,max_mflop)

x = x0;
z = b - A*x;
r_old = precondition(A,z);
rh0 = r_old;
rho_old = 1;
alpha = 1;
om_old = 1;
v_old = x0;
p_old = x0;
t_old = r_old'*r_old;

rn = 1;
bn = 1;
pbn=norm(precondition(A,b));
prn = norm(r_old);

iter = 0;
results=update_results([],'PBICG-STAB',iter,prn);


while (~converged(bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop))

   rho_new = rh0'*r_old;
   ppbeta = (rho_new / rho_old) * (alpha/om_old);
   p_new =  r_old + ppbeta*(p_old - om_old * v_old);

%  z = precondition(A,p_new);
%  v_new = A*z;
   z = A*p_new;
   v_new = precondition(A,z);

   sigma = rh0'*v_new;
   alpha = rho_new / sigma;
   s = r_old - alpha*v_new;

%  z = precondition(A,s);
%  t = A*z;
   z = A*s;
   t = precondition(A,z);

   om_new = t'*s / (t'*t);

   x = x + alpha*p_new + om_new*s;
   r_new = s - om_new*t;

   rho_old = rho_new;
   r_old  = r_new;
   p_old  = p_new;
   v_old  = v_new;
   om_old  = om_new;
   iter = iter+1;
   prn = norm(r_new);
   results=update_results(results,'PBICG-STAB',iter,prn);
end

its=results(:,1);
resids=results(:,4);
