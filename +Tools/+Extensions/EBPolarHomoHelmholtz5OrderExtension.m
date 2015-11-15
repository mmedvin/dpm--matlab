classdef EBPolarHomoHelmholtz5OrderExtension < Tools.Extensions.SuperPolarTwoTupleExtension
 % Equaion Based Homogeneuous Polar Helmholtz 5 Order Extension
    
    
    methods
        function obj = EBPolarHomoHelmholtz5OrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
        end
        
        function res = Expansion(obj,Xi0,Xi1,~)
            [xi0, ~ , xi0tt, ~, xi0tttt   ] = Xi0.Derivatives();
            [xi1, ~ , xi1tt               ] = Xi1.Derivatives();
            
            dr2=obj.dr.^2;
            dr3=obj.dr.^3;
            dr4=obj.dr.^4;
            k2=obj.Coeffs.k.^2;
            
            R2=obj.r0.^2;
            % R3=R.^3;
            R4=obj.r0.^4;
            
            res = (1 - dr2.*k2./2 + dr3.*k2./obj.r0./6  - dr4.*k2.*(3./R2 - k2)/24).*xi0 ...
                +  (obj.dr - dr2./2./obj.r0 + dr3.*(2./R2-k2)./6 - dr4.*(3./R2 - k2)./12./obj.r0).*xi1 ...
                + xi0tt.*dr2.*(-1 + obj.dr./obj.r0 - dr2.*(11./R2 - 2.*k2)./12)./R2./2 ...
                + xi1tt.*dr3.*( obj.dr/obj.r0/2 - 1/3)./R2./2 ...
                + xi0tttt.*dr4./R4./24;

        end
    end
    
end

