%SMOOTH Smooth a vector.
%
%       U_OUT = SMOOTH(LEVEL, B, U, FLAG) applies a smoother defined by the
%       global flag "smooth_flag" and the system AU=B to the vector U on the 
%       given grid level.  FLAG is set to 'pre', 'post', or 'coarse' and
%       defines the number of smoothings applied. 
%
%       Accesses global variables in "include_globals"
%       Accesses global variables in "include_flags"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function   u_out = smooth(level, b, u, flag)

include_globals 
include_flags 

if strcmp(flag, 'pre') == 1
   nu = nu1;
elseif  strcmp(flag, 'post') == 1
   nu = nu2;
elseif strcmp(flag, 'coarse') == 1
   nu = 30;
end

eval(['A = A',num2str(level),';']);

if smooth_flag == WEIGHTED_JACOBI

   D = wt * (1./spdiags(A,[0]));
   for i = 1:nu
       u = u + D.*(b - A*u);
   end

elseif smooth_flag == GAUSS_SEIDEL

   L = tril(A);
   for i = 1:nu
       u = u + L\(b - A*u);
   end

elseif smooth_flag == RB_GAUSS_SEIDEL

   eval(['N = N',num2str(level),';']);   
   red = [1:2:N]; black = [2:2:N];
   D = 1./spdiags(A,[0]);

   for i = 1:nu
      u(red)   = (b(red) - A(red,black) * u(black)) .* D(red);
      u(black) = (b(black) - A(black,red) * u(red)) .* D(black);
   end

end

u_out = u;


