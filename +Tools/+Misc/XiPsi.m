classdef XiPsi < Tools.Common.FunctionWithDerivatives
    
    properties(Access = public)        
        xi0Psi     = 0;
        xi1Psi     = 0;
        xi0PsiTT   = 0;
        xi1PsiTT   = 0;
        xi0PsiTTTT = 0;
    end
    
    methods(Access = public)
        
        function obj = XiPsi(Arguments)
            if nargin > 0
                S = Arguments.Scatterer;
                obj.xi0Psi     = Arguments.PsiBC.xi0Psi(S.th,S.r);
                obj.xi1Psi     = Arguments.PsiBC.xi1Psi(S.th,S.r);
                obj.xi0PsiTT   = Arguments.PsiBC.xi0PsiTT(S.th,S.r);
                obj.xi0PsiTTTT = Arguments.PsiBC.xi0PsiTTTT(S.th,S.r);
                obj.xi1PsiTT   = Arguments.PsiBC.xi1PsiTT(S.th,S.r);
            end
        end
    
        function [xi0,x1,xi0tt,x1tt,x0tttt] = Derivatives(obj)
            
            xi0    = obj.xi0Psi     ;
            x1     = obj.xi1Psi     ;
            xi0tt  = obj.xi0PsiTT   ;
            x1tt   = obj.xi0PsiTTTT ;
            x0tttt = obj.xi1PsiTT   ;

            
        end
        
        
        
    end
end