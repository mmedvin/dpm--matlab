%GET_RHS Generate the right-hand side of a linear system
%
%       B = GET_RHS(NX,NY) generates a vector B of order NX*NY defined
%       by the global flag "rhs_flag".
%
%       Accesses global variables in "include_flags"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function b = get_rhs(nx,ny)

include_flags

N = nx*ny;
[XB,YB] = meshgrid([1:nx]/(nx+1),[1:ny]/(ny+1));
B       = sin(pi*XB).*sin(pi*YB).*sin(sqrt(2)*pi*XB).*sin(sqrt(3)*pi*YB);;
b       = reshape(B,N,1);
