classdef WaveNumberPolarR < Tools.Coeffs.AbstractCoeffs
	% implementation of Polar WaveNumber class k=k(r)
    properties
        k0;
        k;kr;krr;k3r;k4r;k5r;
    end
    
    methods(Static = true)
        function [k,kr] = kkr(r,r0,k0)
            p = (r - r0).*r;
            pr = 2*r-r0;
            if length(r)==1 && r>r0
                p=0;
                pr=0;
            else
                p(r>r0)=0;
                pr(r>r0)=0;
            end
            c=10;
            g = c*p.^6;
            gr = c*6*pr.*(p.^5);
            
            k   = k0.*exp(-g);
            kr  =-k.*gr;
        end
    end
        
    
    methods
        
        function [k,kr,krr,k3r,k4r,k5r] = Derivatives(obj)
            k = obj.k;
            kr=obj.kr;
            krr=obj.krr;
            k3r=obj.k3r;
            k4r=obj.k4r;
            k5r=obj.k5r;
        end
        
        function obj=WaveNumberPolarR(Scatterer,Params)
            try
                r  = Scatterer.R;
            catch
                r = Scatterer.r;
            end
            r0 = Params.r0;
            
           if nargin==2
                p = (r - r0).*r;
                pr = 2*r-r0;
                if length(r)==1 && r>r0
                    p=0;
                    pr=0;
                else
                    p(r>r0)=0;
                    pr(r>r0)=0;
                end
                c=10;
                g = c.*p.^6;
                gr = c.*6.*pr.*(p.^5);
                grr = c.*6.*(5.*(pr.^2) + 2.*p).*(p.^4);
                g3r = c.*60.*(2.*(pr.^3) + 3.*p.*pr).*(p.^3);
                g4r = c.*360.*(pr.^4 + 4.*p.*(pr.^2) + p.^2).*(p.^2);
                g5r = c.*720.*(pr.^5 + 10.*p.*(pr.^3) + 10.*pr.*(p.^2)).*p;
                
                obj.k0=Params.k;
                obj.k   = obj.k0.*exp(-g);
                obj.kr  =-obj.k.*gr;
                obj.krr = obj.k.*(gr.^2 - grr);
                obj.k3r =-obj.k.*(gr.^3 - 3.*gr.*grr + g3r);
                obj.k4r = obj.k.*(gr.^4 - 6.*(gr.^2).*grr + 3.*(grr.^2) + 4.*gr.*g3r - g4r);
                obj.k5r =-obj.k.*(gr.^5 - 10.*(gr.^3).*grr + 10.*(gr.^2).*g3r - 10.*grr.*g3r + 5.*gr.*(3.*(grr.^2) - g4r) + g5r);
            else
                error('Wavenumber constructor called with wrong number of arguments')
            end
        end
        
    end
    
    
end
