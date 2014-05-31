%GET_MATRIX Generate a discrete linear opererator
%
%       [A,N] = GET_MATRIX(NX,NY) generates the matrix A of order N for
%       the problem defined by the global flag "problem_flag" and global
%       parameters "prob_args".
%
%       Accesses global variables in "include_flags"
%       Accesses global variables in "include_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function [A,N] = get_matrix(nx,ny)

include_flags
include_globals

tic

fprintf('Generating %g by %g matrix...',nx*ny,nx*ny)

if (problem_flag == POISSON)

   [A,N] = sp_laplace(nx,ny);

elseif (problem_flag == CUT_SQUARE)

   val=prob_args(1);
   [A,N] = sp_cutsq2d (nx,ny,val);

elseif (problem_flag == HELMHOLTZ)

   lambda=0;
   sigma=prob_args(1);
   [A,N] = sp_convdiff (nx,ny,lambda,sigma);

elseif (problem_flag == POISSON_BOLTZMAN)

   disp(sprintf('POISSON_BOLTZMAN not implemented yet'));

elseif (problem_flag == CONVECT_DIFFUSE)

   lambda=prob_args(1);
   sigma=prob_args(2);
   [A,N] = sp_convdiff (nx,ny,lambda,sigma);

end

toc1=toc;

fprintf(' %g seconds.\n',toc1)
