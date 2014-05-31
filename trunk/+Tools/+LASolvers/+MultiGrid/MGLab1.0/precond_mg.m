%PRECOND_MG Multigrid preconditioner.
%
%       Z = PRECOND_MG(R) applies an MG cycle defined by the global parameters
%       to the vector R.
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function z = precond_mg(r)

level = 1;   
z = mg_cycle(level,r);
