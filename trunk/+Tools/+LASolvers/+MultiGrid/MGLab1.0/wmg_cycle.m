%WMG_CYCLE W-Cycle algorithm.
%
%       U_OUT = WMG_CYCLE(LEVEL, B, U_IN) uses the W-cycle to recursively 
%       solve the linear system AX=B at the given level.  If the optional 
%       starting value U_IN is not passed then U_IN is set to 0's.
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function u_out = wmg_cycle(level, b, u_in)
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
   u_c     = wmg_cycle(level+1, b_c);
   if (level < coarsest(level)),
       u_c     = wmg_cycle(level+1, b_c, u_c);
   end
   correct = interpolate(level, u_c);
   u       = u + correct;
   u_out   = smooth(level, b, u, 'post');
end
