%SP_CUTSQ2D Set up a problem with discontinuous coefficients.
%
%       [A,N] = SP_CUTSQ2D(NX,NY,VAL) sets up the following discretized
%       PDE on an NX by NY grid.  The domain is the unit square and u=0 on 
%       the boundary.
%
%	     - (a*u_x)_x - (a*u_y)_y = f,  
%
%                 / 1 in (0,1) x (0,1) \ (.4,.6) x (.4,.6)
%            a = <
%                 \ val in (.4,.6) x (.4,.6)

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [A,N] = sp_elliptic1(nx,ny,val)
	
N = nx * ny;  N5 = 5 * N; h = 1/(nx+1);
irow = zeros(N5, 1); icol = irow; NZA = irow;

P_XY = zeros(nx+1, ny); 
Q_XY = zeros(nx, ny+1);
S_XY = zeros(nx, ny);

for j = 1:ny
   for k = 1:(nx+1)
      x = (k-0.5)*h; y = j*h;
      if x >= 0.4 & x <= 0.6 & y >= 0.4 & y <= 0.6
         P_XY(k,j) = val;
      else
         P_XY(k,j) = 1;
      end
   end
end

for j = 1:(ny+1)
   for k = 1:nx
      x = k*h; y = (j-0.5)*h;
      if x >= 0.4 & x <= 0.6 & y >= 0.4 & y <= 0.6
         Q_XY(k,j) = val;
      else
         Q_XY(k,j) = 1;
      end
   end
end

for j = 1:ny
   for k = 1:nx
      x = k*h; y = j*h;
      S_XY(k,j) = 0;
   end
end

h2 = h*h;
index = 0;
ii = 0;
for j = 1:ny
   for k = 1:nx
      ii = ii + 1;

      if j ~= 1
         index = index + 1;
         NZA (index) = -Q_XY(k,j);
         irow(index) = ii;
         icol(index) = ii - nx;
      end

      if k ~= 1
         index = index + 1;
         NZA (index) = -P_XY(k,j);
         irow(index) = ii;
         icol(index) = ii - 1;
      end

      index = index + 1;
      NZA (index) = P_XY(k,j) + P_XY(k+1,j) + Q_XY(k,j) + Q_XY(k,j+1) + ...
                    h2 * S_XY(k,j);
      irow(index) = ii;
      icol(index) = ii;

      if k ~= nx
         index = index + 1;
         NZA (index) = -P_XY(k+1,j);
         irow(index) = ii;
         icol(index) = ii + 1;
      end

      if j ~= ny 
         index = index + 1;
         NZA (index) = -Q_XY(k,j+1);
         irow(index) = ii;
         icol(index) = ii + nx;
      end
   end
end            

icol = icol(1:index); irow = irow(1:index);
NZA = NZA(1:index);

A = sparse (irow, icol, NZA, N, N);


      
