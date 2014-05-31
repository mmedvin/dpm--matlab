%SINT   Sine transform.
%       SINT(X) returns the sine transform of vector X.  The length of X
%       must be one less than a power of two.
%       If X is a matrix, the SINT operation is applied to each column.

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function sx = sint(x)

[m, n] = size(x);

y            = zeros(2*(m+1),n);
y(2:m+1,1:n) = x;
z            = fft(y);
sx           = -imag(z(2:m+1,1:n));


