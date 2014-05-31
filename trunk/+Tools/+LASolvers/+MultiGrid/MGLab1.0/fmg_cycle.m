%FMG_CYCLE Full Multigrid Algorithm
%
%       U_OUT = FMG_CYCLE(LEVEL, B) uses the full-multigrid cycle to 
%       recursively solve the linear system AX=B at the given level.
%
%       No global variables are accessed.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function u_out = fmg_cycle(level, b)

if level == coarsest(level)
   u_out   = coarse_grid_solve(level,b);
else 
   b_c	   = restrict(level, b);
   u_c_out = fmg_cycle(level+1, b_c);
   u_f_in  = interpolate(level, u_c_out);
   u_out   = mg_cycle(level, b, u_f_in);
end
