classdef WaveNumberElliptical < Tools.Coeffs.WaveNumberPolarR
    properties
         
        kn;kf;knn;kff;  k3n;k3f;k4n;k4f;knf;knff;knnf;knnff;
        
       % FocalDist;
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
        function [k,kn,kf,knn,kff, k3n,k3f,k4n,k4f,knf,knff,knnf,knnff] = Derivatives(obj)
            k       = obj.k;
            kn      = obj.kn;
            kf      = obj.kf;
            knn     = obj.knn;
            kff     = obj.kff;
            k3n     = obj.k3n;
            k3f     = obj.k3f;
            k4n     = obj.k4n;
            k4f     = obj.k4f;
            knf     = obj.knf;
            knff    = obj.knff;
            knnf    = obj.knnf;
            knnff   = obj.knnff;
        end
        
        function obj = WaveNumberElliptical(Scatterer,Params)%,k0,r0) %(FocalDist,eta,phi,k0,r0)
            % consider reordering of arguments, so it's can be called as constant k....
            obj.IsConstant = false;
            [r,rn,rf,rnn,rff,r3n,r3f,r4n,r4f,rnf,rnff,rnnf,rnnff] = ... 
                Tools.Coeffs.WaveNumberElliptical.chngcoord(Scatterer.FocalDistance,Scatterer.Eta,Scatterer.Phi);
            
            PolarScatterer.r  = r;
            %PolarScatterer.r0 = r0;
            
            obj=obj@Tools.Coeffs.WaveNumberPolarR(PolarScatterer,Params);%r,r0);
            
            obj.kn = obj.kr.*rn ;
            obj.kf = obj.kr.*rf ;
            
            obj.knn = obj.krr.*rn.^2 + obj.kr.*rnn ;
            obj.kff = obj.krr.*rf.^2 + obj.kr.*rff ;
            
            obj.k3n = obj.k3r.*rn.^3 + 3*obj.krr.*rn.*rnn + obj.kr.*r3n;
            obj.k3f = obj.k3r.*rf.^3 + 3*obj.krr.*rf.*rff + obj.kr.*r3f;
            
            obj.k4n = obj.k4r.*rn.^4 + 6*obj.k3r.*rnn.*rn.^2 + 3*obj.krr.*rnn.^2 + 4*obj.krr.*rn.*r3n + obj.kr.*r4n;
            obj.k4f = obj.k4r.*rf.^4 + 6*obj.k3r.*rff.*rf.^2 + 3*obj.krr.*rff.^2 + 4*obj.krr.*rf.*r3f + obj.kr.*r4f;
            
            
            
            obj.knf  = obj.krr.*rf.*rn+obj.kr.*rnf;
            obj.knff = obj.k3r.*rn.*rf.^2 + obj.krr.*rn.*rff + 2*obj.krr.*rf.*rnf + obj.kr.*rnff;
            
            obj.knnf  = obj.k3r.*rf.*rn.^2 + obj.krr.*rf.*rnn + 2*obj.krr.*rn.*rnf + obj.kr.*rnnf;
            obj.knnff = obj.k4r.*(rf.*rn).^2 + obj.k3r.*rff.*rn.^2 + 4*obj.k3r.*rf.*rn.*rnf + obj.k3r.*rnn.*rf.^2 ...
                      + obj.krr.*rff.*rnn + 2*obj.krr.*rf.*rnnf + 2*obj.krr.*rnf.^2 + 2*obj.krr.*rn.*rnff + obj.kr.*rnnff;
            
            
         
            
        
                  

        end
        
        
       
        
        
        
        
        
        
    end
    
    
end
