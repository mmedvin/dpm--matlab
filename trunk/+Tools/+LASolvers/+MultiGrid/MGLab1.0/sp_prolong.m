%SP_PROLONG Generate prolongation matrix
%
%       P = SP_PROLONG(NX1,NY1,NX2,NY2) explicitly generates the bilinear
%       interpolation matrix P for transfering vectors from a coarse (NX2,NY2)
%       grid to a fine (NX1,NY1) grid.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function P = sp_prolong(nx1,ny1,nx2,ny2)

N1=nx1*ny1;   
N2=nx2*ny2;

I=zeros(9*N2,1);   
J=I;   
VAL=I;

v=[1:9]';
for j = 1:ny2
   for i = 1:nx2
 
      k = nx1 + (j-1)*2*nx1 + 2*i; 
      I(v) = [k-nx1-1, k-nx1, k-nx1+1, ...
              k    -1, k,     k    +1, ...
              k+nx1-1, k+nx1, k+nx1+1]';

      k = (j-1)*nx2 + i;      
      J(v) = k * ones(9,1);

      VAL(v) = [1; 2; 1; ...
                2; 4; 2; ...
                1; 2; 1]; 

      v=v+9;

   end
end

P = 0.25 * sparse(I,J,VAL,N1,N2);
