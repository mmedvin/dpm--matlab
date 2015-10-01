classdef EBPolarHelmholtz5OrderExtension < Tools.Extensions.TwoTupleExtension 
 % Equaion Based Polar Helmholtz 5 Order Extension
    
    methods
        function obj = EBPolarHelmholtz5OrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
        end
        
        function res = Expansion(obj,Xi0,Xi1,Src)
            [xi0, xi0t, xi0tt, ~, xi0tttt] = Xi0.Derivatives();
            [xi1, ~   , xi1tt            ] = Xi1.Derivatives();
            
            [f,fr,frr,ftt] = Src.Derivatives();
            [k,kr,krr] = obj.Coeffs.Derivatives('r');
            [~,kt,ktt] = obj.Coeffs.Derivatives('t');
            
            xirr = f - xi1./obj.r0 - xi0tt./(obj.r0.^2) - xi0.*(k.^2);
            xi3r = fr + xi1./(obj.r0.^2) - xirr./obj.r0 + 2*xi0tt./(obj.r0.^3) - xi1tt./(obj.r0.^2) ...
                 - (k.^2).*xi1 - 2.*k.*kr.*xi0;
            
            k2xi0_tt = (k.^2).*xi0tt + 4*k.*kt.*xi0t + (2*(kt.^2)+2*k.*ktt).*xi0 ;
            xittrr = ftt - xi1tt./obj.r0 - xi0tttt./(obj.r0.^2) - k2xi0_tt;%(k.^2).*xi0tt;
            
            xi4r = frr - 2*(kr.^2 + k.*krr).*xi0 - (2./(obj.r0.^3)+ 4*k.*kr).*xi1 ...
                + (2/(obj.r0.^2)-k.^2).*xirr - xi3r./obj.r0 - 6*xi0tt./(obj.r0.^4) ...
                + 4*xi1tt./(obj.r0.^3) - xittrr./(obj.r0.^2);
             
            res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 ;


        end
 
    end
    
end

