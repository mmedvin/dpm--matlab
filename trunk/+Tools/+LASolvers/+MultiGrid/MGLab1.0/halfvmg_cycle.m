%HALFVMG_CYCLE  Half V-cycle Algorithm
%
%       *** NOT IMPLEMENTED YET ***
%
%       U_OUT = HALFVMG_CYCLE(LEVEL, B) uses the half V-Cycle to 
%       recursively solve the linear system AX=B at the given level.
%
%       No global variables are accessed.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function u_out = vmg_cycle(level, b, u_in)

% Use the zero vector for u_in as the default

if nargin == 2,   
   u_in = zeros(size(b));
end

if level == coarsest(level)
   u_out   = coarse_grid_solve(level, b);
else 
   u       = smooth(level, b, u_in, 'pre');
   r       = residual(level, b, u);
   b_c     = restrict(level, r);
   u_c     = vmg_cycle(level+1, b_c);
   correct = interpolate(level, u_c);
   u       = u + correct;
   u_out   = smooth(level, b, u, 'post');
end
