classdef ExLapCrclVarCoeffs < Tools.Exact.SuperExact

    properties(Access = private)
        B; C;
    end
    
    properties%(Access = private)
        u         = 0;
        dudr      = 0; 
        d2udr2    = 0;
        dudtheta  = 0;
        d2udtheta2 = 0;
    end
        
    methods
        
        function obj = ExLapCrclVarCoeffs(Scatterer, Coeffs)
            obj = obj@Tools.Exact.SuperExact(Scatterer, Coeffs);
            
            r0 = Scatterer.r0;
			r = Scatterer.R;
            
            obj.B = Coeffs.B;
            obj.C = Coeffs.C;
            
            p   = r.^2;
            pr  = 2*r;
            prr = 2;
            
            if length(r)==1 && r>r0
                p   = (1 - 1/8/obj.B - 1/obj.B)/4 + ( (r.^4)/2 + r.^2 )/obj.B + obj.C*log(2*r)/obj.B; 

                pr  = 2*( r.^3 + r )/obj.B + obj.C/obj.B./r; 
                prr = 2*( 3*r.^2 + 1 )/obj.B - 2*obj.C/obj.B./r./r;
            else
                p(r>r0)     = (1 - 1/8/obj.B - 1/obj.B)/4 + ( (r(r>r0).^4)/2 + r(r>r0).^2 )/obj.B + obj.C*log(2*r(r>r0))/obj.B; 
                pr(r>r0)    = 2*( r(r>r0).^3 + r(r>r0) )/obj.B + obj.C/obj.B./r(r>r0);
                prr(r>r0)   = 2*( 3*r(r>r0).^2 + 1 )/obj.B - 2*obj.C/obj.B./(r(r>r0).^2);
            end
            obj.u=p;
            obj.dudr=pr;
            obj.d2udr2=prr;
        end
        
        function [u,dudr,d2udr2,dudtheta,d2udtheta2] = Derivatives(obj)
            u         = obj.u;            
            dudr      = obj.dudr;            
            d2udr2    = obj.d2udr2;                    
            dudtheta  = obj.dudtheta;            
            d2udtheta2= obj.d2udtheta2;
        end        
    end
    
end
