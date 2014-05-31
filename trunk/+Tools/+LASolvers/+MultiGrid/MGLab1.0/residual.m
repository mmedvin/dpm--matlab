%RESIDUAL Compute the residual at the given level.
%
%       R = RESIDUAL(LEVEL, B, U) returns the residual R of the system
%       AU=B at the given grid level.
%
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function   r = residual(level, b, u)

include_globals 

eval(['r = b - A', num2str(level), ' * u;']);
