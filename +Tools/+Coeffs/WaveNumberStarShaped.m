classdef WaveNumberStarShaped < Tools.Coeffs.WaveNumberPolarR
    properties
         
        kn;knn;ks;kss;
       FirstTime= 1;
    end 
    
    methods(Static = true)
        function [k,kn] = kkn(FocalDistance,eta,phi,k0,r0)
            
            [r,rn] = Tools.Coeffs.WaveNumberElliptical.chngcoord(FocalDistance,eta,phi);
            
            [k,kr] = Tools.Coeffs.WaveNumberPolarR.kkr(r,r0,k0);
            kn = kr.*rn ;
        end
        
        function [r,rn,rf,rnn,rff,r3n,r3f,r4n,r4f,rnf,rnff,rnnf,rnnff] = chngcoord(FocDist,eta,phi)
                
                x  = FocDist*cosh(eta).*cos(phi);
                y  = FocDist*sinh(eta).*sin(phi);
                xn = FocDist*sinh(eta).*cos(phi);
                yn = FocDist*cosh(eta).*sin(phi);
                xf = -FocDist*cosh(eta).*sin(phi);
                yf = FocDist*sinh(eta).*cos(phi);
                
                
                r = sqrt(x.^2+y.^2);
                
                rn =(xn.*x + yn.*y)./r;
                rf = (xf.*x + yf.*y)./r;
                
                rnn = (cosh(2*eta).*FocDist.^2 - rn.^2)./r;
                rff = (-cos(2*phi).*FocDist.^2  - rf.^2)./r;
                
                r3n =  (4*r.*rn - 3*rn.*rnn)./r;
                r3f = -(4*r.*rf + 3*rf.*rff)./r;
                
                r4n = (4*(r.*rnn + rn.^2 - rn.*r3n) - 3*rnn.^2)./r;
                r4f = (-4*(r.*rff + rf.^2 + rf.*r3f) - 3*rff.^2)./r;
                
                rnf = -rf.*rn./r;
                rnff = -(rff.*rn+2*rf.*rnf)./r;
                rnnf = -(rf.*rnn+2*rn.*rnf)./r;
                rnnff = (-rff.*rnn-2*(rf.*rnnf+rnf.^2+rn.*rnff))./r;
                
                
                %th = atan(tanh(eta).*tan(phi));
                %th = angle(x+1i*y);
                % tn = (tan(phi).*sech(eta).^2 )./ (1 + (tan(phi).*tanh(eta)).^2);
                % tf = (tanh(phi).*sec(eta).^2 )./ (1 + (tan(phi).*tanh(eta)).^2);
                %
                % tnn = -2*sin(2*phi).*sinh(2*eta)./((cos(2*phi)+cosh(2*eta)).^2);
                % tff = -tnn;
                
            end
    end
    
    methods
        function [k,kn,knn,ks,kss,k5r] = Derivatives(obj,Scatterer)
            
            if exist('Scatterer','var')
                if obj.FirstTime
                    obj.FirstTime =0;
                    [x,dx] = Scatterer.XHandle.Derivatives(Scatterer.BasisArg);
                    [y,dy] = Scatterer.YHandle.Derivatives(Scatterer.BasisArg);
                    r_ = sqrt(x.^2 + y.^2);
                    
                    curv = Scatterer.Curvature(Scatterer.BasisArg);
                    [h,ht,htt,h3t] = Scatterer.Metrics(Scatterer.BasisArg);
                    
                    kx = obj.kr.*x./r_;
                    ky = obj.kr.*y./r_;
                    
                    kxx = (obj.kr.*y.^2)./(r_.^3) + (obj.krr.*x.^2)./(r_.^2);
                    kyy = (obj.kr.*x.^2)./(r_.^3) + (obj.krr.*y.^2)./(r_.^2);
                    kxy = -(obj.kr.*x.*y)./(r_.^3) + (obj.krr.*x.*y)./(r_.^2);
                    
                    obj.kn  = (kx.*dy  - ky.*dx)./h;
                    obj.knn = (kxx.*dy.^2 - 2*kxy.*dx.*dy + kyy.*dx.^2 )./(h.^2);
                    
                    obj.ks  = (kx.*dx + ky.*dy)./h;
                    obj.kss = (kxx.*dx.^2 + 2*kxy.*dx.*dy + kyy.*dy.^2)./(h.^2) + curv.*obj.kn;
                end
                
                k       = obj.k;
                kn      = obj.kn;
                knn      = obj.knn;
                ks     = obj.ks;
                kss     = obj.kss;
            else
                [k,kn,knn,ks,kss,k5r] = Derivatives@Tools.Coeffs.WaveNumberPolarR(obj,'r');
                %despite the names it is actually returns [k,kr,krr,k3r,k4r,k5r] here
            end
           
            
           
        end
        
        function obj = WaveNumberStarShaped(Scatterer,Params)          
            obj=obj@Tools.Coeffs.WaveNumberPolarR(Scatterer,Params);
        end
        
        
       
        
        
        
        
        
        
    end
    
    
end
