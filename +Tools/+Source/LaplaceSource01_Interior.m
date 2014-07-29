classdef LaplaceSource01_Interior < Tools.Source.SuperSource
% syms x y c(x,y) d(x,y) a(x,y) b(x,y)
% u(x,y)=sin(c(x,y))*sin(d(x,y));
%    
% diff(a(x,y)*diff(u(x,y),x),x)+diff(b(x,y)*diff(u(x,y),y),y)
% 
% ans =
%  
% diff(a(x, y), x)*(cos(c(x, y))*sin(d(x, y))*diff(c(x, y), x) + cos(d(x, y))*sin(c(x, y))*diff(d(x, y), x)) + diff(b(x, y), y)*(cos(c(x, y))*sin(d(x, y))*diff(c(x, y), y) + cos(d(x, y))*sin(c(x, y))*diff(d(x, y), y)) + a(x, y)*(cos(c(x, y))*sin(d(x, y))*diff(c(x, y), x, x) - sin(c(x, y))*sin(d(x, y))*diff(d(x, y), x)^2 - sin(c(x, y))*sin(d(x, y))*diff(c(x, y), x)^2 + cos(d(x, y))*sin(c(x, y))*diff(d(x, y), x, x) + 2*cos(c(x, y))*cos(d(x, y))*diff(c(x, y), x)*diff(d(x, y), x)) + b(x, y)*(cos(c(x, y))*sin(d(x, y))*diff(c(x, y), y, y) - sin(c(x, y))*sin(d(x, y))*diff(d(x, y), y)^2 - sin(c(x, y))*sin(d(x, y))*diff(c(x, y), y)^2 + cos(d(x, y))*sin(c(x, y))*diff(d(x, y), y, y) + 2*cos(c(x, y))*cos(d(x, y))*diff(c(x, y), y)*diff(d(x, y), y))


%
%syms x y c(y) d(x) a(x,y) b(x,y)
%u(x,y)=sin(c(y)*x)*sin(d(x)*y);
%diff(a(x,y)*diff(u(x,y),x),x)+diff(b(x,y)*diff(u(x,y),y),y)
% 
%ans =
% 
%diff(a(x, y), x)*(y*cos(y*d(x))*sin(x*c(y))*diff(d(x), x) + cos(x*c(y))*sin(y*d(x))*c(y)) + diff(b(x, y), y)*(x*cos(x*c(y))*sin(y*d(x))*diff(c(y), y) + cos(y*d(x))*sin(x*c(y))*d(x)) - a(x, y)*(sin(x*c(y))*sin(y*d(x))*c(y)^2 - y*cos(y*d(x))*sin(x*c(y))*diff(d(x), x, x) + y^2*sin(x*c(y))*sin(y*d(x))*diff(d(x), x)^2 - 2*y*cos(x*c(y))*cos(y*d(x))*c(y)*diff(d(x), x)) - b(x, y)*(sin(x*c(y))*sin(y*d(x))*d(x)^2 - x*cos(x*c(y))*sin(y*d(x))*diff(c(y), y, y) + x^2*sin(x*c(y))*sin(y*d(x))*diff(c(y), y)^2 - 2*x*cos(x*c(y))*cos(y*d(x))*d(x)*diff(c(y), y))



   %  properties (Dependent = true)
   %     Source;%Fn;Ff;Fnn;Fff;    % WN;
   % end
    
    properties(Access = protected)       
        Scatterer; 
		CoeffsClsrHndl;
		CoeffsParams;
		ExParams;
        IsInterior = true;
    end
    
    methods
        function [F,Fn,Ff,Fnn,Fff] = Derivatives(obj)
            
            if obj.IsDummy
                F=0;
                Fn=0;
                Ff=0;
                Fnn=0;
                Fff=0;
			else

                try
                    eta = obj.Scatterer.Eta;
                    phi = obj.Scatterer.Phi;
                    fd  = obj.Scatterer.FocalDistance;
                catch exception
                    if strcmp(exception.identifier,'MATLAB:nonExistentField')
                        eta  = obj.Scatterer.eta;
                        phi = obj.Scatterer.phi;
                        fd  = obj.Scatterer.FocalDistance;
                    else
                        rethrow(exception);
                    end
                end
				
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

				coeffs = obj.CoeffsClsrHndl(obj.Scatterer,obj.CoeffsParams);
				[a,ax,axx,ay,ayy,axy] = coeffs.Derivatives('ax');
                [b,bx,bxx,by,byy,bxy] = coeffs.Derivatives('bx');
                
                if obj.IsInterior
                    [c,cy,cyy,c3y] = obj.ExParams.CoeffsIn.Derivatives('c',x,y);
                    [d,dx,dxx,d3x] = obj.ExParams.CoeffsIn.Derivatives('d',x,y);
                else
                    [c,cy,cyy,c3y] = obj.ExParams.CoeffsOut.Derivatives('c',x,y);
                    [d,dx,dxx,d3x] = obj.ExParams.CoeffsOut.Derivatives('d',x,y);
                end
               
                u   = sin(c.*x).*sin(d.*y);
                ux  = c.*cos(x.*c).*sin(y.*d) + y.*cos(y.*d).*sin(x.*c).*dx;
                uy  = cos(y.*d).*d.*sin(x.*c) + x.*cos(x.*c).*sin(y.*d).*cy;
                uxx = -(c.^2 + (y.*dx).^2).*u + (2.*c.*dx.*cos(x.*c) + sin(x.*c).*dxx).*y.*cos(y.*d);
                uyy = -(d.^2 + (x.*cy).^2).*u + (2.*cy.*d.*cos(y.*d) + sin(y.*d).*cyy).*x.*cos(x.*c);
                
                F = ax.*ux + a.*uxx + by.*uy + b.*uyy;
                
                if nargout>1 
                    
                    uxy = (c.*d + cy.*dx.*x.*y).*cos(c.*x).*cos(d.*y) + dx.*cos(d.*y).*sin(c.*x) + cy.*cos(c.*x).*sin(d.*y) + (-c.*cy.*x - d.*dx.*y).*sin(c.*x).*sin(d.*y);
                    
                    uxyy = -2*c.*cy.*d.*x.*cos(d.*y).*sin(c.*x) + cos(c.*x).*sin(d.*y).*(-2*cy.*d.*dx.*x.*y + cyy) + u.*(-2*d.*dx - 2*(cy.^2).*x - c.*x.*cyy) ...
                           +cos(c.*x).*cos(d.*y).*(2*cy.*d + 2*cy.*dx.*x + dx.*x.*y.*cyy);
                    
                    uyxx = -u.*(2*c.*cy + 2*(dx.^2).*y + d.*dxx.*y) - uy.*(c.^2 + (dx.*y).^2) + (2*c.*dx + 2*cy.*dx.*y + cy.*dxx.*x.*y).*cos(c.*x).*cos(d.*y) ...
                           +(dxx - 2*c.*cy.*dx.*x.*y).*cos(d.*y).*sin(c.*x) - 2*c.*d.*dx.*y.*cos(c.*x).*sin(d.*y);
                    
                    u3x  = -3*dx.*dxx.*u.*(y.^2) + 3*c.*dxx.*y.*cos(c.*x).*cos(d.*y) + (d3x.*y - 2*(c.^2).*dx.*y).*cos(d.*y).*sin(c.*x) ...
                           - 2*c.*((dx.*y).^2).*cos(c.*x).*sin(d.*y);
                    
                    u3y  = -3*cy.*cyy.*u.*(x.^2) + 3*cyy.*d.*x.*cos(c.*x).*cos(d.*y) - 2*(cy.^2).*d.*(x.^2).*cos(d.*y).*sin(c.*x) ...
                           + (c3y.*x - 2*cy.*(d.^2).*x).*cos(c.*x).*sin(d.*y);
                    
                    Fx = a.*u3x + axx.*ux + 2*ax.*uxx + by.*uxy + b.*uxyy + bxy.*uy +   bx.*uyy;

                    Fy = b.*u3y + axy.*ux +   ay.*uxx + a.*uyxx + ax.*uxy + byy.*uy + 2*by.*uyy;
                    
                    %                     m1 = d.^2 + (cy.*x).^2;
                    %
                    %                     Fx = cos(d.*y).*( ...
                    %                         (by.*c.*d + 2*bx.*cy.*d.*x + 4*ax.*c.*dx.*y + 3*a.*c.*dxx.*y + by.*cy.*dx.*x.*y + b.*(2*cy.*(d + dx.*x) + cyy.*dx.*x.*y)).*cos(c.*x) ...
                    %                         + (bxy.*d + 2*ax.*dxx.*y + dx.*(by + axx.*y) + b.*(-2*c.*cy.*d.*x - dx.*(m1).*y) + a.*y.*(d3x - 3*dx.*c.^2 - dx.^3.*(y.^2)) ).*sin(c.*x)) ...
                    %                         + ((axx.*c + bx.*cyy.*x + cy.*(by + bxy.*x) - b.*(-cyy + c.*(m1) + 2*cy.*d.*dx.*x.*y) - a.*(c.^3 + 3*c.*(dx.*y).^2)).*cos(c.*x) ...
                    %                         - (c.*(2*ax.*c + by.*cy.*x) + b.*(2*d.*dx + 2*x.*cy.^2 + c.*cyy.*x) + bx.*(m1) + dx.*y.*(by.*d + 2*ax.*dx.*y + 3*a.*dxx.*y)).*sin(c.*x)).*sin(d.*y);
                    %
                    %                     Fy = sin(c.*x).*((d.*(byy - b.*(d.^2 + 3*(cy.*x).^2)) + ay.*dxx.*y + dx.*(ax + axy.*y) - a.*(d.*c.^2 - dxx + 2*c.*cy.*dx.*x.*y + d.*(dx.*y).^2)).*cos(d.*y) ...
                    %                         - (3*b.*cy.*cyy.*x.^2 + 2*by.*(m1) + ax.*(c.*cy.*x + d.*dx.*y) + a.*(2*c.*cy + 2*y.*dx.^2 + d.*dxx.*y) + ay.*(c.^2 + (dx.*y).^2)).*sin(d.*y))...
                    %                         + cos(c.*x).*(((4*by.*cy + 3*b.*cyy).*d.*x + 2*ay.*c.*dx.*y + ax.*(c.*d + cy.*dx.*x.*y) + a.*(cy.*dxx.*x.*y + 2*dx.*(c + cy.*y))).*cos(d.*y) ...
                    %                         + (axy.*c + ax.*cy + byy.*cy.*x + 2*by.*cyy.*x + b.*x.*(c3y - 3*cy.*d.^2 - cy.^3.*(x.^2)) - a.*(cy.*x.*c.^2  + 2*c.*d.*dx.*y + cy.*x.*(dx.*y).^2)).*sin(d.*y));
                    
                    Fn = Fx.*xn + Fy.*yn;
                    
                end
                
                if nargout>2, Ff = Fx.*xf + Fy.*yf;  end
                
                %if nargout>3, Fnn = 0   ;  end
                
                %if nargout>4, Fff = 0   ;  end
            end
            
        end
            
        function obj = LaplaceSource01_Interior(Scatterer, CoeffsClsrHndl,CoeffsParams,ExParams)
			obj.Scatterer = Scatterer;
			obj.CoeffsClsrHndl = CoeffsClsrHndl;
			obj.CoeffsParams = CoeffsParams;
			obj.ExParams = ExParams;
			obj.IsDummy = false;
		end
        
		function S = Source(obj)
 			S = spalloc(obj.Scatterer.Size(1),obj.Scatterer.Size(2),numel(obj.Scatterer.Np));
			
			%[F,Fn,~,Fnn] = obj.Derivatives();
            [F,Fn] = obj.Derivatives();
			
% 			Coeffs	= obj.CoeffsClsrHndl(obj.Scatterer.TheScatterer,obj.CoeffsParams);
			%Exact	= Tools.Exact.ExLapElps01(obj.Scatterer, obj.ExParams);
			
			
			
			S(obj.Scatterer.Inside) = F(obj.Scatterer.Inside);
            S(1:end,1)=	0; %Exact.u(1:end,1);
            S(1,1:end)= 0; %Exact.u(1,1:end);
            S(1:end,end)= 0; %Exact.u(1:end,end);
            S(end,1:end)= 0; %Exact.u(end,1:end);
			
			S(obj.Scatterer.GridGamma)	= F(obj.Scatterer.GridGamma) ...
										+ obj.Scatterer.deta.*Fn(obj.Scatterer.GridGamma); %...  
										%+ (obj.Scatterer.dr.^2).*Fnn(obj.Scatterer.GridGamma)/2;%taylor
			
		end
	end
end


