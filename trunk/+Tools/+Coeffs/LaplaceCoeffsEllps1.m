classdef LaplaceCoeffsEllps1 < Tools.Coeffs.AbstractCoeffs
	% Implemntation of Laplacian coefficient for DPM
	% in this case sigma=0; a=b;
	properties
		B;
		
        ax;axx;ay;ayy;axy;
        bx;bxx;by;byy;bxy;
        
		a; an; ann; %a3n; a4n; a5n;
           af; aff;
		b; bf; bff;
           bn; bnn;
           bnf;% bff; %b3f; b4f; b5f;
		%br;btr;
		sigma; sigma_r; sigma_rr; sigma_3r; sigma_4r; sigma_5r;
	end
	
	methods
		
		function [d,dr,drr,d3r,d4r,d5r] = Derivatives(obj,WhichOne)
            %d=0;dr=0;drr=0;d3r=0;d4r=0;d5r=0;
			switch WhichOne
				case {'a','an'}
					d   = obj.a;
					dr  = obj.an;
					drr = obj.ann;
					%d3r = obj.a3r;
					%d4r = obj.a4r;
					%d5r = obj.a5r;
				case 'af'
					d   = obj.a;
					dr  = obj.af;
					drr = obj.aff;
					%d3r = obj.a3r;
					%d4r = obj.a4r;
					%d5r = obj.a5r;                    
                case 'ax'
                    d   = obj.a;
                    dr  = obj.ax;
                    drr = obj.axx;
                    d3r = obj.ay;
                    d4r = obj.ayy;
                    d5r = obj.axy;
				case 'bn'
					d   = obj.b;
					dr  = obj.bn;
					drr = obj.bnn;                    
					%d3r = obj.b3t;
					%d4r = obj.b4t;
					%d5r = obj.b5t;
				case {'b','bf'}
					d   = obj.b;
					dr  = obj.bf;
					drr = obj.bff;                    
					%d3r = obj.b3t;
					%d4r = obj.b4t;
					%d5r = obj.b5t;
                    
                case 'bx'
                    d   = obj.b;
                    dr  = obj.bx;
                    drr = obj.bxx;
                    d3r = obj.by;
                    d4r = obj.byy;
                    d5r = obj.bxy;
				case 'sigma'
					d   = obj.sigma;
					dr  = obj.sigma_r;
					drr = obj.sigma_rr;
					d3r = obj.sigma_3r;
					d4r = obj.sigma_4r;
					d5r = obj.sigma_5r;
				otherwise
					error('unknown coeff')
			end
		end
		
		function obj=LaplaceCoeffsEllps1(Scatterer,Params)
			
 %           try
                fd   = Scatterer.FocalDistance;
                eta  = Scatterer.Eta;
                phi  = Scatterer.Phi;
                %eta0 = Scatterer.Eta0;
%            catch
%                 fd = Params.FocalDistance;
%                 
%                 r  = Scatterer.r;
%                 th = Scatterer.th;
%                 
%                 ez = acosh( r.*(cos(th) + 1i*sin(th))/fd );
%                 eta = real(ez);
%                 phi = imag(ez);
%            end
            
            c = Params.ca;
            d = Params.da;
            e = Params.ea;
				
            x  = fd*cosh(eta).*cos(phi);
            y  = fd*sinh(eta).*sin(phi);
            xn = fd*sinh(eta).*cos(phi); %xnn = x;
            yn = fd*cosh(eta).*sin(phi); %ynn=y;   
            xf =-fd*cosh(eta).*sin(phi); %xff=-x;
            yf = fd*sinh(eta).*cos(phi); %yff=-y;
            
			obj.a   = e + sin(c*x+d*y)/2;
            
            obj.ax = c*cos(c*x+d*y)/2;
            obj.ay = d*cos(c*x+d*y)/2;
            obj.axx=-(c^2)*sin(c*x+d*y)/2;
            obj.ayy=-(d^2)*sin(c*x+d*y)/2;
            obj.axy=-c*d*sin(c*x+d*y)/2;
            
            obj.an  = obj.ax.*xn + obj.ay.*yn;
			obj.ann = obj.axx.*xn.^2 + obj.ayy.*yn.^2 + 2*obj.axy.*xn.*yn + obj.ax.*x + obj.ay.*y;
            
            obj.af  = obj.ax.*xf + obj.ay.*yf;
            obj.aff = obj.axx.*(xf.^2) + obj.ayy.*(yf.^2) + 2*obj.axy.*xf.*yf + obj.ax.*(-x) + obj.ay.*(-y);
            
			%obj.a3n = 0;  obj.a4n = 0;  obj.a5n = 0;
            
            if Params.WithB                
                c              = Params.cb;
                d              = Params.db;
                e              = Params.eb;
                
                obj.b   = e + cos(c*x+d*y)/2;
                
                obj.bx = -c*sin(c*x+d*y)/2;
                obj.by = -d*sin(c*x+d*y)/2;
                obj.bxx=-(c^2)*cos(c*x+d*y)/2;
                obj.byy=-(d^2)*cos(c*x+d*y)/2;
                obj.bxy=-c*d*cos(c*x+d*y)/2;
            else
                obj.b=obj.a;
                obj.bx=obj.ax;
                obj.bxx=obj.axx;
                obj.by=obj.ay;
                obj.byy=obj.ayy;
                obj.bxy=obj.axy;
            end

			obj.bn  = obj.bx.*xn + obj.by.*yn;
            obj.bnn = obj.bxx.*xn.^2 + obj.byy.*yn.^2 + 2*obj.bxy.*xn.*yn + obj.bx.*x + obj.by.*y;
            
            obj.bf  = obj.bx.*xf + obj.by.*yf;
            obj.bff = obj.bxx.*xf.^2 + obj.byy.*yf.^2 + 2*obj.bxy.*xf.*yf + obj.bx.*(-x) + obj.by.*(-y);
           
			obj.bnf = obj.bxx.*xn.*xf + obj.byy.*yn.*yf + obj.bxy.*(xn.*yf + xf.*yn) + obj.bx.*(-y) + obj.by.*x;
            
			obj.sigma=0; obj.sigma_r=0; obj.sigma_rr=0; obj.sigma_3r=0; obj.sigma_4r=0; obj.sigma_5r=0;
		end
	end
end
