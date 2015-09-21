classdef EBPolarHomoHelmholtz7OrderExtension < Tools.Extensions.SuperPolarTwoTupleExtension
    % Equaion Based Homogeneuous Polar Helmholtz 7 Order Extension
    
    
    methods
        function obj = EBPolarHomoHelmholtz7OrderExtension(Arguments)
            obj = obj@Tools.Extensions.SuperPolarTwoTupleExtension(Arguments);
        end
        
        function res = Expansion(obj,Xi0,Xi1,~)
            
            [xi0, ~, xi0tt, ~, xi0tttt, xi0tttttt   ] = Xi0.Derivatives();
            [xi1, ~, xi1tt, ~, xi1tttt              ] = Xi1.Derivatives();
            
            [k,kr,krr,k3r] = obj.Coeffs.Derivatives('r');
 
            f=0; fr=0; frr=0; f3r=0; ftt=0; frtt=0; ftttt=0;fttrr=0;
            
            xirr = f - xi1./obj.r0 - xi0tt./(obj.r0.^2) - xi0.*(k.^2);
            xi3r = fr + xi1./(obj.r0.^2) - xirr./obj.r0 + 2.*xi0tt./(obj.r0.^3) - xi1tt./(obj.r0.^2) - (k.^2).*xi1 - 2.*k.*kr.*xi0;
            xittrr = ftt - xi1tt./obj.r0 - xi0tttt./(obj.r0.^2) - (k.^2).*xi0tt;
            xi4r = frr - 2.*(kr.^2 + k.*krr).*xi0 - (2./(obj.r0.^3)+ 4.*k.*kr).*xi1 + (2./(obj.r0.^2)-k.^2).*xirr - xi3r./obj.r0 - 6.*xi0tt./(obj.r0.^4) + 4.*xi1tt./(obj.r0.^3) - xittrr./(obj.r0.^2);
            
            xitt3r = frtt + xi1tt./(obj.r0.^2) - xittrr./obj.r0 + 2.*xi0tttt./(obj.r0.^3) - xi1tttt./(obj.r0.^2)-...
                (k.^2).*xi1tt - 2.*k.*kr.*xi0tt;
            
            xi5r = f3r - 2.*(3.*kr.*krr+k.*k3r).*xi0 - 6.*(-1./(obj.r0./4)+kr.^2+k.*krr).*xi1-...
                6.*(1./(obj.r0.^3) + k.*kr).*xirr + (3./(obj.r0.^2) - k.^2).*xi3r - xi4r./obj.r0 +...
                24.*xi0tt./(obj.r0.^5) - 18.*xi1tt./(obj.r0.^4)+ 6.*xittrr./(obj.r0.^3)  - xitt3r./(obj.r0.^2);
            
            %                res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 + (obj.dr.^5).*xi5r/120;
            
            xi4trr = ftttt - xi1tttt./obj.r0 - xi0tttttt./(obj.r0.^2) - xi0tttt.*(k.^2);
            xitt4r = fttrr - 2.*xi1tt./(obj.r0.^3) + (2./(obj.r0.^2)-k.^2).*xittrr - xitt3r./obj.r0 - 6*xi0tttt./(obj.r0.^4) + 4*xi1tttt./(obj.r0.^3) - xi4trr./(obj.r0.^2);
            
            xi6r = 24*xirr./(obj.r0.^4) - 24*xi1./(obj.r0.^5) - 12*xi3r./(obj.r0.^3) + (4./(obj.r0.^2) - k^2).*xi4r - xi5r./obj.r0 ...
                - 120*xi0tt./(obj.r0.^6) + 96*xi1tt./(obj.r0.^5) - 36*xittrr./(obj.r0.^4) + 8*xitt3r./(obj.r0.^3) - xitt4r./(obj.r0.^2);
            
            res = xi0 + obj.dr.*xi1 + (obj.dr.^2).*xirr/2 + (obj.dr.^3).*xi3r/6 + (obj.dr.^4).*xi4r/24 + (obj.dr.^5).*xi5r/120 + (obj.dr.^6).*xi6r/720;
            
        end
    end
    
end

