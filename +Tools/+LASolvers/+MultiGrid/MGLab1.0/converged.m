%CONVERGED  Return whether the method has converged
%
%       C = CONVERGED (BN,PBN,RN,PRN,ITER,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP)
%       Performs a series of tests to determine whether the calling iterative
%       method has converged.  
%
%         BN and PBN are the 2-norms of the original system and the 
%            preconditioned system respectively.
%         RN and PRN are the 2-norms of the residuals of the system and the 
%            preconditioned system respectively.
%         ITER is the current iteration number.
%         RTOL and PRTOL are the tolerances for the residuals RN and PRN
%            respectively.
%         MAX_IT is the limit on the number of iterations.
%         MAX_TIME is the limit on the time (not implemented yet).
%         MAX_MFLOP is the limit on the number of megaflops (not implemented 
%            yet).
%         

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 7 June 1995

function c = converged (bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop)

    c = (( rtol~=0 &  rn <= rtol*bn) ...
       | (prtol~=0 & prn <= prtol*pbn) ...
       | (iter >= max_it));
