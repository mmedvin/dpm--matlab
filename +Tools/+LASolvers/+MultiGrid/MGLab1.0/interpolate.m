%INTERPOLATE Transfer correction from the coarse grid to the current grid.
%
%       CORRECT = INTERPOLATE(LEVEL, U_C_OUT) uses the interpolation scheme
%       defined by "interp_flag" to transfer the vector U_C_OUT on the coarse 
%       level LEVEL+1 to the vector CORRECT on the current level LEVEL.
%
%       Accesses global variables in "include_flags"
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function correct = interpolate(level, u_c_out)

include_flags
include_globals 
extract_globals

IX_c = 2:(nx_c+1); IY_c = 2:(ny_c+1);
IX_f = 2:(nx_f+1); IY_f = 2:(ny_f+1);

if interp_flag == LINEAR

   COARSE(IX_c,IY_c) = reshape(u_c_out,nx_c,ny_c);
   FINE = interp2(X_c,Y_c,COARSE,X_f,Y_f);
   correct = reshape(FINE(IX_f,IY_f), N_f, 1);

%  COARSE = reshape(u_c_out,nx_c,ny_c); 
%  COARSE = [zeros(1,nx_c+2); ...  
%            zeros(ny_c,1) COARSE zeros(ny_c,1); ...
%            zeros(1,nx_c+2)];
%  FINE = zeros(nx_f,ny_f);
%  IX  = 2:2:nx_f-1;
%  IX2 = 1:2:nx_f;
%  IY  = 2:2:ny_f-1;
%  IY2 = 1:2:ny_f;
%  ix = 2:nx_c+1; ixm=[1 ix];  ixp=[ix nx_c+2];
%  iy = 2:ny_c+1; iym=[1 iy];  iyp=[iy ny_c+2];
%  FINE(IX,IY) = COARSE(ix,iy);
%  FINE(IX,IY2) = 0.5*(COARSE(ix,iym)+COARSE(ix,iyp));
%  FINE(IX2,IY) = 0.5*(COARSE(ixm,iy)+COARSE(ixp,iy));
%  FINE(IX2,IY2) = 0.25*(COARSE(ixm,iym)+COARSE(ixp,iym)...
%                       +COARSE(ixm,iyp)+COARSE(ixp,iyp));
 
%  correct = reshape(FINE, N_f, 1);

elseif interp_flag == CUBIC

   COARSE(IX_c,IY_c) = reshape(u_c_out,nx_c,ny_c);
   FINE = interp2(X_c,Y_c,COARSE,X_f,Y_f,'cubic');
   correct = reshape(FINE(IX_f,IY_f), N_f, 1);

elseif interp_flag == EXPLICIT_BILINEAR

   eval(['PROLONG = ARRAY',num2str(level), ';']);
   correct = PROLONG * u_c_out;

end
