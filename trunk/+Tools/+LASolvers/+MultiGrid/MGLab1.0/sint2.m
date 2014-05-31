%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function W = sint2(U)
%
%	See fft2
%	.' does a non-Hermitian transpose
%
W = sint(sint(U.').');
