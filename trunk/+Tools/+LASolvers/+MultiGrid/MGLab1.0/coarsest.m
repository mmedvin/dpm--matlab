%COARSEST Return the number of the coarsest grid level
%
%       COARSEST(level) returns the global variable "coarse_level".  The
%       levels are numbered such that the finest mesh is 1.
%
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function clevel = coarsest(level)

include_globals

clevel = coarse_level;
