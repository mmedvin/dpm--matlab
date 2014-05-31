%SP_CONVDIFF Set up a convection-diffusion problem.
%
%       [A,N] = SP_CONVDIFF(NX,NY,LAM,SIGMA) sets up the following discretized
%       PDE on an NX by NY grid.  The domain is the unit square and u=0 on 
%       the boundary.
%
%	     - u_xx - u_yy + lam u_x + sigma u = f
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [A,N] = sp_convdiff(nx,ny,lam,sigma)
N = nx * ny;  N5 = 5 * N;
irow = zeros(N5, 1); icol = irow; NZA = irow;

h = 1/(nx+1);  h2 = h*h;

index = 0;
ii = 0;
for j = 1:ny
   for k = 1:nx
      ii = ii + 1;

      if j ~= 1
         index = index + 1;
         NZA (index) = -1.0;
         irow(index) = ii;
         icol(index) = ii - nx;   % S
      end

      if k ~= 1
         index = index + 1;
         NZA (index) = -1.0 - lam * h;
         irow(index) = ii;
         icol(index) = ii - 1;    % W
      end

      index = index + 1;
      NZA (index) = 4.0  + lam * h + sigma * h2;
      irow(index) = ii;
      icol(index) = ii;           % P (self)

      if k ~= nx
         index = index + 1;
         NZA (index) = -1.0;
         irow(index) = ii;
         icol(index) = ii + 1;    % E
      end

      if j ~= ny 
         index = index + 1;
         NZA (index) = -1.0;
         irow(index) = ii;
         icol(index) = ii + nx;   % N
      end
   end
end            

icol = icol(1:index); irow = irow(1:index);
NZA = NZA(1:index);

A = sparse (irow, icol, NZA, N, N);


      

