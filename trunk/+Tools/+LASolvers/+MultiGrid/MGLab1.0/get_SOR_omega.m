%GET_SOR_OMEGA Return the SOR relaxation parameter from globals
%
%       OMEGA = GET_SOR_OMEGA returns the SOR relaxation parameter.
%
%       Accesses global variables in "include_globals"
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 20 April 1995

function om = get_SOR_omega

include_globals

om = SOR_omega;
