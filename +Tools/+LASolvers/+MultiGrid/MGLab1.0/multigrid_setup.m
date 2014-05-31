%MULTIGRID_SETUP Generate the linear systems for the coarse multigrid levels
%
%       MULTIGRID_SETUP generates the coarse-grid linear systems for the 
%       currently defined problem.
%
%       Accesses global variables in "include_flags"
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function multigrid_setup

include_flags
include_globals

[X1,Y1] = meshgrid([0:nx1+1]/(nx1+1),[0:ny1+1]/(ny1+1));
ARRAY1 = zeros(nx1+2, ny1+2); 

if coarse_level >= 2
   nx2 = (nx1+1)/2 - 1;
   ny2 = (ny1+1)/2 - 1;
   [A2,N2] = get_matrix(nx2,ny2);
   [X2,Y2] = meshgrid([0:nx2+1]/(nx2+1),[0:ny2+1]/(ny2+1));
   if interp_flag == EXPLICIT_BILINEAR
      ARRAY1 = sp_prolong(nx1,ny1,nx2,ny2);
   else
      ARRAY2 = zeros(nx2+2, ny2+2);
   end
end

if coarse_level >= 3
   nx3 = (nx2+1)/2 - 1;
   ny3 = (ny2+1)/2 - 1;
   [A3,N3] = get_matrix(nx3,nx3);
   [X3,Y3] = meshgrid([0:nx3+1]/(nx3+1),[0:ny3+1]/(ny3+1));
   if interp_flag == EXPLICIT_BILINEAR
      ARRAY2 = sp_prolong(nx2,nx2,nx3,nx3);
   else
      ARRAY3 = zeros(nx3+2, ny3+2);
   end
end

if coarse_level >= 4
   nx4 = (nx3+1)/2 - 1;
   ny4 = (ny3+1)/2 - 1;
   [A4,N4] = get_matrix(nx4,nx4);
   [X4,Y4] = meshgrid([0:nx4+1]/(nx4+1),[0:ny4+1]/(ny4+1));
   if interp_flag == EXPLICIT_BILINEAR
      ARRAY3 = sp_prolong(nx3,nx3,nx4,nx4);
   else
      ARRAY4 = zeros(nx4+2,ny4+2);
   end
end

if coarse_level >= 5
   nx5 = (nx4+1)/2 - 1;
   ny5 = (ny4+1)/2 - 1;
   [A5,N5] = get_matrix(nx5,nx5);
   [X5,Y5] = meshgrid([0:nx5+1]/(nx5+1),[0:ny5+1]/(ny5+1));
   if interp_flag == EXPLICIT_BILINEAR
      ARRAY4 = sp_prolong(nx4,nx4,nx5,nx5); 
   else
      ARRAY5 = zeros(nx5+2,ny5+2);
   end
end
