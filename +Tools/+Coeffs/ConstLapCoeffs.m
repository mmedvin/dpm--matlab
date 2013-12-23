classdef ConstLapCoeffs < Tools.Coeffs.AbstractCoeffs
    properties
        a; ar;
        b; br;
        sigma;
    end
    
%     methods(Static = true)
%         function [k,kr] = kkr(r,r0,k0)
%             k=obj.k;
%             kr=0;
%         end
%     end
        
    
    methods
        
        function [a,ar,b,br,sigma] = Derivatives(obj)
            a = obj.a;
            ar = obj.ar;
            b = obj.b;  
            br = obj.br;
            sigma = obj.sigma;
        end
        
        function obj=ConstLapCoeffs(~,Params) %(k0,r,r0)
            obj.a       = Params.a;
            obj.ar=0;
            obj.b       = Params.b;
            obj.br=0;
            obj.sigma   = Params.sigma;
        end
        
    end
    
    
end
