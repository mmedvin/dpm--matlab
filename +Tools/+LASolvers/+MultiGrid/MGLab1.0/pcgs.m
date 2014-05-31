%PCGS   Preconditioned conjugate gradient squared method
%
%       [X,RESIDS,ITS]=PCGS(A,B,X0,RTOL,PRTOL,MAX_IT,MAX_TIME,MAX_MFLOP)
%       solves the system AX = B using the preconditioned conjugate gradient 
%       squared method with the given tolerances and limits.
%

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 8 June 1995

function [XNEW,resids,its] = pcgs (A,b,x0,rtol,prtol,max_it,max_time,max_mflop);

    XOLD = x0;

%%% initialize iteration

    Z = b - A*XOLD;
    ROLD  = precondition(A,Z);
    
    QOLD = zeros(size(x0));    POLD = QOLD;
    RHOOLD = 1;
    R0 = ROLD;

%%% initialize iters, resids


    bn = 1;
    rn = 1;
    pbn = precondition(A,b);
    prn = norm(ROLD);

    iter    = 0;
    results = update_results([],'CGS',iter,prn);

%%% begin loop 


     while (~converged(bn,pbn,rn,prn,iter,rtol,prtol,max_it,max_time,max_mflop))

       RHONEW = R0' * ROLD;                                 % inner
       beta = RHONEW / RHOOLD;
       UNEW = ROLD + beta*QOLD;                             % saxpy
       PNEW = UNEW + beta*(QOLD + beta*POLD);               % saxpy,saxpy

       Z = A * PNEW;       VNEW = precondition(A,Z);        % inv(M)*A

       sigma = R0' * VNEW;       % (see footnote (!))       % inner
       alpha = RHONEW  / sigma; 
       QNEW = UNEW - alpha * VNEW;                          % saxpy
       VNEW = alpha * (UNEW + QNEW);                        % saxpy
       XNEW = XOLD + VNEW;                                  % vecadd

       Z = A * VNEW;       Z = precondition(A,Z);           % inv(M)*A

       RNEW = ROLD - Z;                                     % vecadd

       %%% prepare for next iteration 

       RHOOLD = RHONEW;       ROLD  = RNEW;
       QOLD  = QNEW;          POLD  = PNEW;
       XOLD  = XNEW;          prn = norm(RNEW);

       iter = iter + 1;
       results=update_results(results,'CGS',iter,prn);

    end

its=results(:,1);
resids=results(:,4);


%%% -- footnote ----------------------------
%%% ! SIAM J. Sci. Stat. Comput.  Vol.10 No.1 p.44 1989 says 
%%%   R0' * VNEW should be R0' * UNEW.  I think this is a typing mistake
%%%   (SIAM J. Sci. Stat. Comput.  Vol.13 No.2 p.632 1992 says
%%%    R0' * VNEW though)
