%MAX_LEVEL Return the largest number of levels given the grid size.
%
%       LEVEL = MAX_LEVEL(NX) returns the maximum number of grid levels 
%       LEVEL for a grid with NX as the number of grid points along the 
%       shortest axis.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function level = max_level(nx)

level = log2(nx+1);

