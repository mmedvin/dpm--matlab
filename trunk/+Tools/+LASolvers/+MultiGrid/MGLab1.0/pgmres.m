%PCGMRES Preconditioned generalized minimum residual (GMRES) method.
%
%       [X,RESIDS,ITS]=PCGMRES(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP,...
%       RESTART) solves the system AX = B using the restarted preconditioned 
%       GMRES method with the given tolerances, limits and restart parameter.
%
 
% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [x,resids,its] = pcgmres(A,b,x0,...
     rtol,prtol,max_it,max_time,max_mflop,restart)

   x_old = x0;

   r_old = b - A*x_old;
   r_old = precondition(A,r_old);
   prn = norm(r_old);

   B = zeros(restart+1,1);

%%% begin loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   iter = 0;
   reit     = restart+1;
   results=update_results([],'PGMRES',iter,prn);

   bn = 1;
   rn = 1;
   pbn=norm(precondition(A,b));

   while (~converged(bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop))

      %%% check for restart

      if (reit > restart)
          clear V H; 

          reit = 1;
          X0   = x_old;
          B(1) = prn;
          V(:,1) = r_old / prn;
      end

      %%% basic iteration 

      v_temp = A*V(:,reit);
      v_temp = precondition(A,v_temp);

      for j = 1:reit
         H(j,reit) = V(:,j)' * v_temp;
         v_temp = v_temp - (H(j,reit) * V(:,j));
      end

      H(reit+1,reit) = sqrt(v_temp' * v_temp);
      V(:,reit+1) = v_temp / H(reit+1,reit);
   
      Y = H(1:reit+1,1:reit) \ B(1:reit+1);

      little_res = B(1:reit+1) - H(1:reit+1,1:reit) * Y(1:reit);
      prn = norm(little_res);

      x_new = X0 + V(:,1:reit) * Y(1:reit);
      r_old = V(:,1:reit+1) * little_res;

      x_old  = x_new;

      iter = iter + 1;
      reit = reit + 1;
      results=update_results(results,'PGMRES',iter,prn);
      
   end
 
its=results(:,1);
resids=results(:,4);

x=x_new;
