classdef ExLapElps01 < Tools.Exact.SuperExact

    properties(Access = private)
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
            error(' ')
			x = FocalDist*cosh(eta).*cos(phi);
			y = FocalDist*sinh(eta).*sin(phi);
			
			e = x.^2 - y.^2;
			
			if length(eta)==1 && eta > Eta0
				e = sin(c*x).*sin(d*y);
			else
				e(eta > Eta0) = sin(x(eta > Eta0)).*cos(y(eta > Eta0));
			end
		end			
	
	function dnde = dndExact(FocalDist,eta,phi,Eta0)
        error(' ')
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
        error(' ')
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
        
        function ICu = InterfaceConditionU(obj,phi)
            
            x  = obj.FocalDistance*cosh(obj.Eta0).*cos(phi);
            y  = obj.FocalDistance*sinh(obj.Eta0).*sin(phi);
            xn = obj.FocalDistance*sinh(obj.Eta0).*cos(phi);
            yn = obj.FocalDistance*cosh(obj.Eta0).*sin(phi);
            
            
            c1 = obj.Coeffs.CoeffsIn.Derivatives('c',x,y);
            d1 = obj.Coeffs.CoeffsIn.Derivatives('d',x,y);

            c0 = obj.Coeffs.CoeffsOut.Derivatives('c',x,y);
            d0 = obj.Coeffs.CoeffsOut.Derivatives('d',x,y);            
            
            ICu = sin(c1.*x).*sin(d1.*y) - sin(c0.*x).*sin(d0.*y);
                        
        end
        
        function ICun = InterfaceConditionUn(obj,phi)
                              
            x  = obj.FocalDistance*cosh(obj.Eta0).*cos(phi);
            y  = obj.FocalDistance*sinh(obj.Eta0).*sin(phi);
            xn = obj.FocalDistance*sinh(obj.Eta0).*cos(phi);
            yn = obj.FocalDistance*cosh(obj.Eta0).*sin(phi);
            
            
            [c1,c1y] = obj.Coeffs.CoeffsIn.Derivatives('c',x,y);
            [d1,d1x] = obj.Coeffs.CoeffsIn.Derivatives('d',x,y);

            [c0,c0y] = obj.Coeffs.CoeffsOut.Derivatives('c',x,y);
            [d0,d0x] = obj.Coeffs.CoeffsOut.Derivatives('d',x,y);            
                        

            u0x = y.*d0x.*sin(c0.*x).*cos(d0.*y) + c0.*cos(c0.*x).*sin(d0.*y);
            u0y = x.*c0y.*cos(c0.*x).*sin(d0.*y) + d0.*cos(d0.*y).*sin(c0.*x);
            
            u1x = y.*d1x.*sin(c1.*x).*cos(d1.*y) + c1.*cos(c1.*x).*sin(d1.*y);
            u1y = x.*c1y.*cos(c1.*x).*sin(d1.*y) + d1.*cos(d1.*y).*sin(c1.*x);
                        
            ICun = (u1x.*xn + u1y.*yn) - (u0x.*xn + u0y.*yn);
            
        end
        
        function [ICu,ICun] = InterfaceCondition(obj,N)
            
            phi=linspace(0,2*pi, N+1);
            phi = phi(1:N);
            
            x  = obj.FocalDistance*cosh(obj.Eta0).*cos(phi);
            y  = obj.FocalDistance*sinh(obj.Eta0).*sin(phi);
            xn = obj.FocalDistance*sinh(obj.Eta0).*cos(phi);
            yn = obj.FocalDistance*cosh(obj.Eta0).*sin(phi);
            
            
            [c1,c1y] = obj.Coeffs.CoeffsIn.Derivatives('c',x,y);
            [d1,d1x] = obj.Coeffs.CoeffsIn.Derivatives('d',x,y);

            [c0,c0y] = obj.Coeffs.CoeffsOut.Derivatives('c',x,y);
            [d0,d0x] = obj.Coeffs.CoeffsOut.Derivatives('d',x,y);            

            
            ICu = sin(c1.*x).*sin(d1.*y) - sin(c0.*x).*sin(d0.*y);
            

            u0x = y.*d0x.*sin(c0.*x).*cos(d0.*y) + c0.*cos(c0.*x).*sin(d0.*y);
            u0y = x.*c0y.*cos(c0.*x).*sin(d0.*y) + d0.*cos(d0.*y).*sin(c0.*x);
            
            u1x = y.*d1x.*sin(c1.*x).*cos(d1.*y) + c1.*cos(c1.*x).*sin(d1.*y);
            u1y = x.*c1y.*cos(c1.*x).*sin(d1.*y) + d1.*cos(d1.*y).*sin(c1.*x);
                        
            ICun = (u1x.*xn + u1y.*yn) - (u0x.*xn + u0y.*yn);
            
        end
        
        function [u,ux,uxx,uy,uyy,uxy] = exact(obj,x,y,c,cy,cyy,d,dx,dxx)
            u   = sin(c.*x).*sin(d.*y);
            ux  = c.*cos(x.*c).*sin(y.*d) + y.*cos(y.*d).*sin(x.*c).*dx;
            uy  = cos(y.*d).*d.*sin(x.*c) + x.*cos(x.*c).*sin(y.*d).*cy;
            uxx = -(c.^2 + (y.*dx).^2).*u + (2.*c.*dx.*cos(x.*c) + sin(x.*c).*dxx).*y.*cos(y.*d);
            uyy = -(d.^2 + (x.*cy).^2).*u + (2.*cy.*d.*cos(y.*d) + sin(y.*d).*cyy).*x.*cos(x.*c);
            uxy = (c.*d + cy.*dx.*x.*y).*cos(c.*x).*cos(d.*y) + dx.*cos(d.*y).*sin(c.*x) + cy.*cos(c.*x).*sin(d.*y) + (-c.*cy.*x - d.*dx.*y).*sin(c.*x).*sin(d.*y);
        end
        
        function obj = ExLapElps01(Scatterer, Coeffs)		
            obj = obj@Tools.Exact.SuperExact(Scatterer, Coeffs);
            
			obj.FocalDistance   = Scatterer.FocalDistance;
            obj.Eta             = Scatterer.Eta;
            obj.Phi             = Scatterer.Phi;			
            obj.Eta0            = Scatterer.Eta0;
           
            fd = obj.FocalDistance;
            eta = obj.Eta;
            phi = obj.Phi;
            
            x  = fd*cosh(eta).*cos(phi);
            y  = fd*sinh(eta).*sin(phi);
            xn = fd*sinh(eta).*cos(phi);
            %xnn = x;
            yn = fd*cosh(eta).*sin(phi);
            %ynn=y;
            xf =-fd*cosh(eta).*sin(phi);
            %xff=-x;
            yf = fd*sinh(eta).*cos(phi);
            %yff=-y;           
            
            [c,cy,cyy] = Coeffs.CoeffsIn.Derivatives('c',x,y);
            [d,dx,dxx] = Coeffs.CoeffsIn.Derivatives('d',x,y);
            

            [obj.u,ux,uxx,uy,uyy,uxy] = exact(obj,x,y,c,cy,cyy,d,dx,dxx);
            
            obj.dudeta = ux.*xn + uy.*yn;
            %obj.dudphi = ux.*xf + uy.*yf;
            
            obj.d2udeta2 = uxx.*xn.^2 + uyy.*yn.^2 + 2*uxy.*xn.*yn;
            %obj.d2udphi2 = uxx.*xf.^2 + uyy.*yf.^2 + uxy.*xf.*yf;
            
            if length(eta)==1 && eta > obj.Eta0
                [c,cy,cyy] = Coeffs.CoeffsOut.Derivatives('c',x,y);
                [d,dx,dxx] = Coeffs.CoeffsOut.Derivatives('d',x,y);
                
                [obj.u,ux,uxx,uy,uyy,uxy] = exact(obj,x,y,c,cy,cyy,d,dx,dxx);
                
                obj.dudeta = ux.*xn + uy.*yn;
                %obj.dudphi = ux.*xf + uy.*yf;
                
                obj.d2udeta2 = uxx.*xn.^2 + uyy.*yn.^2 + 2*uxy.*xn.*yn;
                %obj.d2udphi2 = uxx.*xf.^2 + uyy.*yf.^2 + uxy.*xf.*yf;
            else
                [c,cy,cyy] = Coeffs.CoeffsOut.Derivatives('c',x,y);
                [d,dx,dxx] = Coeffs.CoeffsOut.Derivatives('d',x,y);
                
                if size(c)==1
                    c = c*ones(size(eta));
                    cy = cy*ones(size(eta));
                    cyy = cyy*ones(size(eta));
                end

                if size(d)==1
                    d   = d*ones(size(eta));
                    dx  = dx*ones(size(eta));
                    dxx = dxx*ones(size(eta));
                end
                                             
                c   = c(eta > obj.Eta0); 
                cy  = cy(eta > obj.Eta0);
                cyy = cyy(eta > obj.Eta0);
                
                d   = d(eta > obj.Eta0);
                dx  = dx(eta > obj.Eta0);
                dxx = dxx(eta > obj.Eta0);
                
                x=x(eta > obj.Eta0);
                y=y(eta > obj.Eta0);
                xn=xn(eta > obj.Eta0);
                yn=yn(eta > obj.Eta0);
                
                [obj.u(eta > obj.Eta0),ux,uxx,uy,uyy,uxy] = exact(obj,x,y,c,cy,cyy,d,dx,dxx);
                
                obj.dudeta(eta > obj.Eta0) = ux.*xn + uy.*yn;
                %obj.dudphi = ux.*xf + uy.*yf;
                
                obj.d2udeta2(eta > obj.Eta0) = uxx.*xn.^2 + uyy.*yn.^2 + 2*uxy.*xn.*yn;
                %obj.d2udphi2 = uxx.*xf.^2 + uyy.*yf.^2 + uxy.*xf.*yf;
            end

        end
        
        function [u,dudeta,d2udeta2,dudphi,d2udphi2] = Derivatives(obj)
            u         = obj.u;            
            dudeta    = obj.dudeta;      
            d2udeta2  = obj.d2udeta2;
			if nargout>3
				error('higher derivative not implemented')
			end
        end        
    end
    
end
