%============================
% Multigrid V-Cycle Algorithm
%============================

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function u_out = demo2_Vcycle(level, b, u_in, u_true, nx,iter)

% Use the zero vector for u_in as the default

   if level == 1 & iter == 1
      subplot(1,1,1), hold off, cla
      e = u_true - u_in; E = reshape(e, nx, nx);
      subplot(2,2,1), hold off, surf(E);  hold on
       title(['Initial error'])
      subplot(2,2,2), hold off, surf(abs(sint2(E))); hold on
       title('Absolute value of initial error in Fourier space')
      pause(4)
      %%%%%%%%%%mypause, print -deps D10.eps
   end

if nargin == 2,   
   u_in = zeros(size(b));
end

if level == coarsest(level)
   u_out   = coarse_grid_solve(level, b);
else 
   u       = smooth(level, b, u_in, 'pre');
   r       = residual(level, b, u);
   b_c     = restrict(level, r);
   u_c     = demo2_Vcycle(level+1, b_c, zeros(size(b_c)), ...
                            u_true, nx,iter);
   correct = interpolate(level, u_c);
   u       = u + correct;
   if level == 1
      e = u_true - u; E = reshape(e, nx, nx);
      subplot(2,2,1), hold off, surf(E); hold on, 
       title(['Error after coarse grid correction, iter = ', num2str(iter)])
      subplot(2,2,2), hold off, surf(abs(sint2(E))); hold on
       title('Absolute value of error in Fourier space')
      pause(3)
   end
   u_out   = smooth(level, b, u, 'post');
   if level == 1
      e = u_true - u_out; E = reshape(e, nx, nx);
      subplot(2,2,3), hold off, surf(E); hold on
       title(['Error after post-smoothing, iter = ', num2str(iter)])
      subplot(2,2,4), hold off, surf(abs(sint2(E))); hold on
       title('Absolute value of error in Fourier space')
      pause(3)
      %%%%%%%%%%mypause, eval(['print -deps D1',num2str(iter),'.eps'])
   end
end


