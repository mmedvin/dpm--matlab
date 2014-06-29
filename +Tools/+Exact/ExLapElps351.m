classdef ExLapElps351 < Tools.Exact.SuperExact

    properties(Access = private)
        B; C;
		
		FocalDistance;
        Eta;
        Phi;
		Eta0;

    end
    
    properties%(Access = private)
        u           = 0;
        dudeta      = 0; 
        d2udeta2    = 0;
%         dudphi      = 0;
%         d2udphi2    = 0;
	end
        
	methods(Static = true)
		
		function e = Exact(FocalDist,eta,phi,Eta0)
			x = FocalDist*cosh(eta).*cos(phi);
			y = FocalDist*sinh(eta).*sin(phi);
			
			e = x.^2 - y.^2;
			
			if length(eta)==1 && eta > Eta0
				e = sin(x).*cos(y);
			else
				e(eta > Eta0) = sin(x(eta > Eta0)).*cos(y(eta > Eta0));
			end
		end			
	
	function dnde = dndExact(FocalDist,eta,phi,Eta0)
		x  = FocalDist*cosh(eta).*cos(phi);
		y  = FocalDist*sinh(eta).*sin(phi);
		xn = FocalDist*sinh(eta).*cos(phi);
		yn = FocalDist*cosh(eta).*sin(phi);
		
		dnde = 2*x.*xn - 2.*y.*yn;
		
		if length(eta)==1 && eta > Eta0
			dnde = cos(x).*cos(y).*xn -  sin(x).*sin(y).*yn;
		else
			dnde(eta > Eta0) = cos(x(eta > Eta0)).*cos(y(eta > Eta0)).*xn(eta > Eta0) -  sin(x(eta > Eta0)).*sin(y(eta > Eta0)).*yn(eta > Eta0);
		end
		
		
	end
    
    function d2nde2 = d2ndExact2(FocalDist,eta,phi,Eta0)
        x  = FocalDist*cosh(eta).*cos(phi);
        y  = FocalDist*sinh(eta).*sin(phi);
        xn = FocalDist*sinh(eta).*cos(phi);
        yn = FocalDist*cosh(eta).*sin(phi);
        xnn=x;
        ynn=y;
        
        d2nde2 = 2*xn.^2 + 2*x.^2 - (2*yn.^2 + 2.*y.^2);
        
        if length(eta)==1 && eta > Eta0
            d2nde2 = cos(x).*cos(y).*xnn - sin(x).*cos(y).*(xn.^2 +  yn.^2 + 2*xn.*yn)   -  sin(x).*sin(y).*ynn;
        else
            d2nde2(eta > Eta0)  = cos(x(eta > Eta0)).*cos(y(eta > Eta0)).*xnn(eta > Eta0) ...
                                - sin(x(eta > Eta0)).*cos(y(eta > Eta0)).*(xn(eta > Eta0).^2 +  yn(eta > Eta0).^2 + 2*xn(eta > Eta0).*yn(eta > Eta0)) ...
                                - sin(x(eta > Eta0)).*sin(y(eta > Eta0)).*ynn(eta > Eta0);
        end
        
        
    end
    
    
end
	
    methods
        
        function obj = ExLapElps351(Scatterer, Coeffs)		
            obj = obj@Tools.Exact.SuperExact(Scatterer, Coeffs);
            
			obj.FocalDistance   = Scatterer.FocalDistance;
            obj.Eta         = Scatterer.Eta;
            obj.Phi         = Scatterer.Phi;			
            obj.Eta0 = Scatterer.Eta0;
			
			obj.u  = obj.Exact(obj.FocalDistance,obj.Eta,obj.Phi,obj.Eta0);
			obj.dudeta = obj.dndExact(obj.FocalDistance,obj.Eta,obj.Phi,obj.Eta0);
            obj.d2udeta2 = obj.d2ndExact2(obj.FocalDistance,obj.Eta,obj.Phi,obj.Eta0);
        end
        
        function [u,dudeta,d2udeta2,dudphi,d2udphi2] = Derivatives(obj)
            u         = obj.u;            
            dudeta    = obj.dudeta;      
			if nargout>0
				error('higher derivative not implemented')
			end
        end        
    end
    
end
