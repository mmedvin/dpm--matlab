function P = PsiBC(Param,Exact,Fn,t)
    Ft = 1;
    if nargin == 4
        Ft = Fn(t);
    end
    
            P.xi0Psi     =@(theta,r) Exact.Psi             (theta,Param(r)).*Ft;
            P.xi1Psi     =@(theta,r) Exact.DrDPsi          (theta,Param(r)).*Ft;
            P.xi0PsiTT   =@(theta,r) Exact.DPsiDThetaTheta (theta,Param(r)).*Ft;
            P.xi1PsiTT   =@(theta,r) Exact.DPsiDrThetaTheta(theta,Param(r)).*Ft;
            P.xi0PsiTTTT =@(theta,r) Exact.DPsiD4Theta     (theta,Param(r)).*Ft;
end
