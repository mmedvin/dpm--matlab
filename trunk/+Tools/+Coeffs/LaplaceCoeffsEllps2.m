classdef LaplaceCoeffsEllps2 < Tools.Coeffs.AbstractCoeffs
	% Implemntation of Laplacian coefficient for DPM
	% in this case sigma=0; a=b;
	properties
		B;
		
        ax;axx; a3x; ay;ayy; a3y; axy; axyy; axxy;
        bx;bxx;b3x; by;byy;b3y; bxy; bxxy; bxyy;
        
		a; an; ann; a3n; 
        af; aff; 
        anf; anff;
        
		b; bf; bff; b3f; 
        bn; bnn;
        bfn; bfnn;
		
		sigma; sigma_r; sigma_rr; sigma_3r; sigma_4r; sigma_5r;
	end
	
	methods
        
		function [d,dr,drr,d3r,d4r,d5r,d6r,d7r] = Derivatives(obj,WhichOne)
            %d=0;dr=0;drr=0;d3r=0;d4r=0;d5r=0;
			switch WhichOne
				case 'a'
					d   = obj.a;
					dr  = obj.an;
					drr = obj.ann;
					d3r = obj.a3n;
					d4r = obj.af;
					d5r = obj.aff;
                    d6r = obj.anf;
                    d7r = obj.anff;
                case 'ax'
                    d   = obj.a;
                    dr  = obj.ax;
                    drr = obj.axx;
                    d3r = obj.ay;
                    d4r = obj.ayy;
                    d5r = obj.axy;
				case 'b'
					d   = obj.b;
					dr  = obj.bf;
					drr = obj.bff;                    
					d3r = obj.b3f;
					d4r = obj.bn;
                    d5r = obj.bnn;
					d6r = obj.bfn;
                    d7r = obj.bfnn;                   
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
		
		function obj=LaplaceCoeffsEllps2(Scatterer,Params)
			
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
            xn = fd*sinh(eta).*cos(phi); 
            xnf=-fd*sinh(eta).*sin(phi); 
            xnn = x;
            x3n=xn;
            yn = fd*cosh(eta).*sin(phi); 
            ynf = fd*cosh(eta).*cos(phi); 
            ynn=y;           
            y3n=yn;
            xf =-fd*cosh(eta).*sin(phi); 
            xnnf = xf;
            xff=-x;
            x3f=-xf;
            xffn=-xn;
            yf = fd*sinh(eta).*cos(phi); 
            ynnf=yf;
            yff=-y;
            y3f=-yf;
            yffn=-ynf;
            
			obj.a   = c + d*exp(e*x);
            
            obj.ax = d*e*exp(e*x);
            obj.ay = 0;
            obj.axx= d*(e^2)*exp(e*x);
            obj.ayy= 0;
            obj.axy= 0;
            obj.a3x= d*(e^3)*exp(e*x);
            obj.a3y= 0;
            obj.axyy= 0;
            obj.axxy= 0;
            
            
            %[a,an,ann, a3n,af,aff, anf,anff]
            obj.an  = obj.ax.*xn + obj.ay.*yn;
			obj.ann = obj.axx.*(xn.^2) + obj.ayy.*(yn.^2) + 2*obj.axy.*xn.*yn + obj.ax.*xnn + obj.ay.*ynn;
            obj.a3n = obj.a3x.*(xn.^3) + obj.a3y.*(yn.^3) + 3*obj.axyy.*xn.*(yn.^2) + 3*obj.axxy.*(xn.^2).*yn + 2*obj.axx.*xn.*xnn + 2*obj.ayy.*yn.*ynn ...
                    + 3*obj.axy.*(xnn.*yn + xn.*ynn) + obj.ax.*x3n + obj.ay.*y3n;
            
            obj.af  = obj.ax.*xf + obj.ay.*yf;
            obj.aff = obj.axx.*xf.^2 + obj.ayy.*yf.^2 + 2*obj.axy.*xf.*yf + obj.ax.*xff + obj.ay.*yff;
			
            obj.anf = obj.axx.*xn.*xf + obj.ayy.*yn.*yf + obj.axy.*(xn.*yf + xf.*yn) + obj.ax.*xnf + obj.ay.*ynf;            
            obj.anff= obj.a3y.*(yf.^2).*yn + obj.axyy.*(xn.*(yf.^2) + 2*xf.*yf.*yn) + obj.axxy.*(2*xf.*yf.*xn + (xf.^2).*yn) + obj.a3x.*(xf.^2).*xn  ...
                    + obj.axx.*(xn.*xff+2*xf.*xnf) +  obj.ayy.*(yn.*yff + 2*yf.*ynf) + obj.axy.*(2*xnf.*yf + 2*xf.*ynf +yn.*xff + xn.*yff) ...
                    + obj.ax.*xffn + obj.ay.*yffn;
            
            if Params.WithB                
                c              = Params.cb;
                d              = Params.db;
                e              = Params.eb;
                
                obj.b   = c + d*exp(e*y);
                
                obj.bx = 0;
                obj.by = d*e*exp(e*y);
                obj.bxx= 0;
                obj.byy= d*(e^2)*exp(e*y);
                obj.bxy= 0;
                obj.b3x= 0;
                obj.b3y= d*(e^3)*exp(e*y);
                obj.bxyy= 0;
                obj.bxxy= 0;
            else
                error('unexpected initialization')
                obj.b=obj.a;
                obj.bx=obj.ax;
                obj.bxx=obj.axx;
                obj.by=obj.ay;
                obj.byy=obj.ayy;
                obj.bxy=obj.axy;
            end
            %obj.a=obj.b;
            %obj.ax=obj.bx;
            %obj.axx=obj.bxx;
            %obj.ay=obj.by;
            %obj.ayy=obj.byy;
            %obj.axy=obj.bxy;
            %obj.an  = obj.ax.*xn + obj.ay.*yn;
            %obj.ann = obj.axx.*xn.^2 + obj.ayy.*yn.^2 + 2*obj.axy.*xn.*yn;

            %[b,bf,bff,b3f,bn,bnn,bfn,bfnn] 
			obj.bf	= obj.bx.*xf + obj.by.*yf;
            obj.bff = obj.bxx.*(xf.^2) + obj.byy.*(yf.^2) + 2*obj.bxy.*xf.*yf  + obj.bx.*xff + obj.by.*yff;
            obj.b3f = obj.b3x.*(xf.^3) + obj.b3y.*(yf.^3) + obj.bxyy.*3*xf.*(yf.^2) + obj.bxxy.*3*(xf.^2).*yf + 2*obj.bxx.*xf.*xff + 2*obj.byy.*yf.*yff ...
                    + 3*obj.bxy.*(xff.*yf + xf.*yff) + obj.bx.*x3f + obj.by.*y3f;

			obj.bn  = obj.bx.*xn + obj.by.*yn;
            obj.bnn = obj.bxx.*xn.^2 + obj.byy.*yn.^2 + 2*obj.bxy.*xn.*yn + obj.bx.*xnn + obj.by.*ynn;
            
			obj.bfn = obj.bxx.*xn.*xf + obj.byy.*yn.*yf + obj.bxy.*(xn.*yf + xf.*yn) + obj.bx.*xnf + obj.by.*ynf;
            
			obj.bfnn= obj.b3x.*(xn.^2).*xf + obj.b3y.*(yn.^2).*yf + obj.bxyy.*(xf.*(yn.^2) + 2*yf.*xn.*yn) + obj.bxxy.*(2*xf.*xn.*yn + (xn.^2).*yf) ...
                    + obj.bxx.*(xf.*xnn+2*xn.*xnf) +  obj.byy.*(yf.*ynn + 2*yn.*ynf) + obj.bxy.*(2*xnf.*yn + 2*xn.*ynf +yf.*xnn + xf.*ynn) ...
                    + obj.bx.*xnnf + obj.by.*ynnf;
            
			obj.sigma=0; obj.sigma_r=0; obj.sigma_rr=0; obj.sigma_3r=0; obj.sigma_4r=0; obj.sigma_5r=0;
		end
	end
end
