%RESTRICT Transfer residual from the current grid to the next coarser grid.
%
%       RHS_C = RESTRICT(LEVEL,R) uses the restriction scheme defined by 
%       "restrict_flag" to transfer the vector R on the current level LEVEL 
%       to the vector RHS_C on the next coarser level LEVEL+1.
%
%       Accesses global variables in "include_flags"
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function rhs_c = restrict(level,r)

include_globals 
extract_globals
include_flags 

%  2-D RESTRICTIONS

   nx0_f = nx_f+2;
   ny0_f = ny_f+2;
   N0_f = nx0_f*ny0_f;
   dx=1;
   dy=nx0_f;

%  Generate r0 by padding r with boundary elements (0's)

   r0 = zeros(N0_f,1);
   for iy=1:ny_f
   for ix=1:nx_f
       r0(nx0_f+1 + ix + nx0_f*(iy-1)) = r(ix+nx_f*(iy-1));
   end
   end

%  Generate indicies of corresponding coarse vector elements in fine vector

   I = zeros(N_c,1);
   for iy=1:ny_c
   for ix=1:nx_c
       I(ix + nx_c*(iy-1)) = 2*ix + 2*iy*nx0_f + 1;
   end
   end

if restrict_flag == INJECTION

   rhs_c = r0(I);

elseif restrict_flag == HALF_WEIGHTING

   rhs_c = .5*r0(I) + ...
          .125*(r0(I+dx) + r0(I-dx) + r0(I+dy) + r0(I-dy));

elseif restrict_flag == FULL_WEIGHTING

   rhs_c = .25*r0(I) + ...
          .125*(r0(I+dx) + r0(I-dx) + r0(I+dy) + r0(I-dy)) + ...
          .0625*(r0(I+dx+dy) + r0(I-dx+dy) + r0(I+dx-dy) + r0(I-dx-dy));

elseif restrict_flag == BILINEAR_ADJOINT

   eval(['PROLONG = ARRAY',num2str(level), ';']);
   rhs_c = PROLONG' * r;

end

rhs_c = 4*rhs_c;

