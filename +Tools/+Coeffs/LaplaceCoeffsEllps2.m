classdef LaplaceCoeffsEllps2 < Tools.Coeffs.AbstractCoeffs
	% Implemntation of Laplacian coefficient for DPM
	% in this case sigma=0; a=b;
	properties
		B;
		
        ax;axx; a3x; ay;ayy; a3y; axy; axyy; axxy;
        bx;bxx;b3x; by;byy;b3y; bxy; bxxy; bxyy;
        
		a; an; ann; a3n; 
            af; aff; a3f;
        anf; anff; annf;
        
		b; bf; bff; b3f; 
           bn; bnn; b3n;
           bnf; bnnf; bnff;
		
		sigma; sigma_r; sigma_rr; sigma_3r; sigma_4r; sigma_5r;
	end
	
	methods
        
		function [d,dr,drr,d3r,d4r,d5r,d6r,d7r,d8r,d9r] = Derivatives(obj,WhichOne)
            %d=0;dr=0;drr=0;d3r=0;d4r=0;d5r=0;
			switch WhichOne
				case 'a'
					d   = obj.a;
					dr  = obj.an;
					drr = obj.ann;
					d3r = obj.a3n;
					d4r = obj.af;
					d5r = obj.aff;
                    d6r = obj.a3f;
                    d7r = obj.anf;
                    d8r = obj.anff;                    
                    d9r = obj.annf;

                case 'ax'
                    d   = obj.a;
                    dr  = obj.ax;
                    drr = obj.axx;
                    d3r = obj.ay;
                    d4r = obj.ayy;
                    d5r = obj.axy;
				case 'b'
					d   = obj.b;
					dr  = obj.bn;
					drr = obj.bnn;
					d3r = obj.b3n;
					d4r = obj.bf;
					d5r = obj.bff;
                    d6r = obj.b3f;
                    d7r = obj.bnf;
                    d8r = obj.bnff;                    
                    d9r = obj.bnnf;
     
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
            yn = fd*cosh(eta).*sin(phi);
            xf =-yn;%fd*cosh(eta).*sin(phi); 
            yf = xn;%fd*sinh(eta).*cos(phi);
            
            ynf = x;
            xnf=-y; 
            xnn = x;
            x3n=xn;
            ynn=y;           
            y3n=yn;
            xnnf = xf;
            xff=-x;
            x3f=-xf;
            xffn=-xn;
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
			obj.a3f = obj.a3x.*(xf.^3) + obj.a3y.*(yf.^3) + 3*obj.axyy.*xf.*(yf.^2) + 3*obj.axxy.*(xf.^2).*yf + 2*obj.axx.*xf.*xff + 2*obj.ayy.*yf.*yff ...
                    + 3*obj.axy.*(xff.*yf + xf.*yff) + obj.ax.*x3f + obj.ay.*y3f;
            
			
            obj.anf = obj.axx.*xn.*xf + obj.ayy.*yn.*yf + obj.axy.*(xn.*yf + xf.*yn) + obj.ax.*xnf + obj.ay.*ynf;            
            obj.anff= obj.a3y.*(yf.^2).*yn + obj.axyy.*(xn.*(yf.^2) + 2*xf.*yf.*yn) + obj.axxy.*(2*xf.*yf.*xn + (xf.^2).*yn) + obj.a3x.*(xf.^2).*xn  ...
                    + obj.axx.*(xn.*xff+2*xf.*xnf) +  obj.ayy.*(yn.*yff + 2*yf.*ynf) + obj.axy.*(2*xnf.*yf + 2*xf.*ynf +yn.*xff + xn.*yff) ...
                    + obj.ax.*xffn + obj.ay.*yffn;
            obj.annf= (obj.a3x.*xf + obj.axxy.*yf).*(xn.^2) + 2*xn.*xnf.*obj.axx ... 
                    + (obj.axyy.*xf + obj.a3y.*yf).*(yn.^2) + 2*yn.*ynf.*obj.ayy ...
                    + (obj.axxy.*xf + obj.axyy.*yf).*xn.*yn + 2*obj.axy.*(xnn.*yn + xn.*ynn) ...
                    + (obj.axx.*xf + obj.axy.*yf).*xnn + obj.ax.*xnnf ...
                    + (obj.axy.*xf + obj.ayy.*yf).*ynn + obj.ay.*ynnf ;
    
            
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
            obj.b3n = obj.b3x.*(xn.^3) + obj.b3y.*(yn.^3) + 3*obj.bxyy.*xn.*(yn.^2) + 3*obj.bxxy.*(xn.^2).*yn + 2*obj.bxx.*xn.*xnn + 2*obj.byy.*yn.*ynn ...
                    + 3*obj.bxy.*(xnn.*yn + xn.*ynn) + obj.bx.*x3n + obj.by.*y3n;
            
            
            
			obj.bnf = obj.bxx.*xn.*xf + obj.byy.*yn.*yf + obj.bxy.*(xn.*yf + xf.*yn) + obj.bx.*xnf + obj.by.*ynf;
            
            obj.bnff= obj.b3y.*(yf.^2).*yn + obj.bxyy.*(xn.*(yf.^2) + 2*xf.*yf.*yn) + obj.bxxy.*(2*xf.*yf.*xn + (xf.^2).*yn) + obj.b3x.*(xf.^2).*xn  ...
                    + obj.bxx.*(xn.*xff+2*xf.*xnf) +  obj.byy.*(yn.*yff + 2*yf.*ynf) + obj.bxy.*(2*xnf.*yf + 2*xf.*ynf +yn.*xff + xn.*yff) ...
                    + obj.bx.*xffn + obj.by.*yffn;
                
            obj.bnnf= (obj.b3x.*xf  + obj.bxxy.*yf).*(xn.^2)+ 2*xn.*xnf.*obj.bxx ... 
                    + (obj.bxyy.*xf + obj.b3y.*yf).*(yn.^2) + 2*yn.*ynf.*obj.byy ...
                    + (obj.bxxy.*xf + obj.bxyy.*yf).*xn.*yn + 2*obj.bxy.*(xnn.*yn + xn.*ynn) ...
                    + (obj.bxx.*xf  + obj.bxy.*yf).*xnn     + obj.bx.*xnnf ...
                    + (obj.bxy.*xf  + obj.byy.*yf).*ynn     + obj.by.*ynnf ;
            
			obj.sigma=0; obj.sigma_r=0; obj.sigma_rr=0; obj.sigma_3r=0; obj.sigma_4r=0; obj.sigma_5r=0;
		end
	end
end
