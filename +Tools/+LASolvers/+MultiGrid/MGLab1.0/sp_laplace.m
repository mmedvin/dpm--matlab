%SP_LAPLACE Set up the discrete Laplace equation.
%
%       [A,N] = SP_LAPLACE(NX,NY) sets up the following discretized
%       PDE on an NX by NY grid.  The domain is the unit square and u=0 on 
%       the boundary.
%
%	     - u_xx - u_yy = f
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [A,N] = sp_laplace(nx,ny)

N = nx * ny;  
N5 = 5 * N;
irow = zeros(N5, 1);
icol = irow; 
NZA = irow;

index = 0;
row = 0;

for j = 1:ny
   for k = 1:nx

      row = row + 1;

      if j > 1
         index = index + 1;
         NZA (index) = -1.0;
         irow(index) = row;
         icol(index) = row - nx;   % South
      end

      if k > 1
         index = index + 1;
         NZA (index) = -1.0;
         irow(index) = row;
         icol(index) = row - 1;    % West
      end

      index = index + 1;
      NZA (index) = 4.0;
      irow(index) = row;
      icol(index) = row;           % P (self)

      if k < nx
         index = index + 1;
         NZA (index) = -1.0;
         irow(index) = row;
         icol(index) = row + 1;    % East
      end

      if j < ny 
         index = index + 1;
         NZA (index) = -1.0;
         irow(index) = row;
         icol(index) = row + nx;   % North
      end
   end
end            

icol = icol(1:index); 
irow = irow(1:index);
NZA = NZA(1:index);

A = sparse (irow, icol, NZA, N, N);
