%MG_CYCLE Multigrid cycle algorithm
%
%       U_OUT = MG_CYCLE(LEVEL, B, U_IN) uses the multigrid cycle defined
%       by the global variable "cycle_flag" to recursively solve the linear 
%       system AX=B at the given level.  If the optional starting value U_IN 
%       is not passed then U_IN is set to 0's.
%
%       Accesses global variables in "include_flags"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function u_out = mg_cycle(level, b, u_in)

include_flags

% Use the zero vector for u_in as the default

if nargin == 2,   
   u_in = zeros(size(b));
end

if (cycle_flag == V_CYCLE)
    u_out = vmg_cycle(level, b, u_in);
elseif (cycle_flag == W_CYCLE)
    u_out = wmg_cycle(level, b, u_in);
elseif (cycle_flag == HALF_V_CYCLE)
    u_out = halfvmg_cycle(level, b, u_in);
end
