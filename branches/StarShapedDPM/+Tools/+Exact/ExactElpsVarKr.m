classdef ExactElpsVarKr < Tools.Exact.SuperExact

    properties%(Access = private)
        u           = 0;
        dudeta      = 0; 
%         d2udeta2    = 0;
%         dudphi      = 0;
%         d2udphi2    = 0;
    end
    
    
    properties(Access = private)
%         Scatterer;
        FocalDistance;
        Eta;
        Phi;
%         WaveNumber;
    end
    
    methods(Static = true)
        
        function e = Exact(FocalDist,eta,phi,k)
            x = FocalDist*cosh(eta).*cos(phi);
            e = exp(1i.*k.*x);
        end
        
        function dnde = dndExact(FocalDist,eta,phi,k,kn)
            e = Tools.Exact.ExactElpsVarKr.Exact(FocalDist,eta,phi,k);
            dnde = 1i.*FocalDist.*cos(phi).*(kn.*cosh(eta)+k.*sinh(eta)).*e;
        end
    end    
    
    methods
        
        function obj = ExactElpsVarKr(Scatterer, WaveNumber)
            obj = obj@Tools.Exact.SuperExact(Scatterer, WaveNumber);
            
            obj.FocalDistance   = Scatterer.FocalDistance;
            obj.Eta         = Scatterer.Eta;
            obj.Phi         = Scatterer.Phi;
            
            obj.u  = obj.Exact(obj.FocalDistance,obj.Eta,obj.Phi,obj.WaveNumber.k);
            obj.dudeta = obj.dndExact(obj.FocalDistance,obj.Eta,obj.Phi,obj.WaveNumber.k,obj.WaveNumber.kn);                       
        end
        
        function [u,dudeta,d2udeta2,dudphi,d2udphi2] = Derivatives(obj)
            u           = obj.u;
            
            dudeta      = obj.dudeta;
            
            d2udeta2    = 1i*obj.FocalDistance*obj.u.*cos(obj.Phi) ...
                        .*(obj.WaveNumber.knn.*cosh(obj.Eta)   ...
                        + 2.*obj.WaveNumber.kn.*sinh(obj.Eta)  ...
                        + obj.WaveNumber.k.*cosh(obj.Eta))     ...
                        + 1i*obj.FocalDistance*obj.dudeta.*cos(obj.Phi)...
                        .*(obj.WaveNumber.kn.*cosh(obj.Eta) + obj.WaveNumber.k.*sinh(obj.Eta));
                    
            dudphi      = 1i*obj.FocalDistance*obj.u .* cosh(obj.Eta).*(obj.WaveNumber.kf.*cos(obj.Phi) - obj.WaveNumber.k.*sin(obj.Phi));
            
            d2udphi2    = 1i*obj.FocalDistance * obj.u .* cosh(obj.Eta).*(obj.WaveNumber.kff .* cos(obj.Phi) ...
                        - 2*obj.WaveNumber.kf.*sin(obj.Phi)-obj.WaveNumber.k.*cos(obj.Phi)) ...
                        + 1i*obj.FocalDistance * dudphi .* cosh(obj.Eta).*(obj.WaveNumber.kf .* cos(obj.Phi) ...
                        - obj.WaveNumber.k .* sin(obj.Phi));
        end
        
        function [s,sn,snn,sf,sff] = calc_s(obj)
            
%             [k,kn,kf,knn,kff,k3n,k3f,k4n,k4f,knf,knff,knnf,knnff] =
%             elptc_k(n,f,k0);
            
            Metrics =Tools.Metrics.EllipticalMetrics(obj.FocalDistance,obj.Eta,obj.Phi);          
            [h,hn,hnn,~,~,hf,hff] = Metrics.metrics();
            [k,kn,kf,knn,kff, k3n,k3f,k4n,k4f,knf,knff,knnf,knnff] = obj.WaveNumber.Derivatives();
            
            
            FD = obj.FocalDistance;
            n = obj.Eta;
            f = obj.Phi;
            
            
            s   = (1i*FD*(-2*kf.*cosh(n).*sin(f) +  2*kn.*sinh(n).*cos(f) + (kff + knn).*cosh(n).*cos(f))...
                - (FD^2)*(((kn.*cosh(n) + k.*sinh(n)).*cos(f)).^2 + ((kf.*cos(f) - k.*sin(f)).*cosh(n)).^2))./(h.^2) + k.^2;
            
            if nargout > 1
                sn  = ( 2*(k.^2).* (-(FD^2)*cosh(n).*sinh(n) + h.*hn) + k.*((FD^2)*sin(2*f).*sinh(2*n).*kf - 2*((-h.^2 + (FD^2)*(cosh(2*n) ...
                    + cos(2*f).*sinh(n).^2)).*kn + (FD^2)*cos(f).*cosh(n).*(-cosh(n).*sin(f).*knf + cos(f).*sinh(n).*knn))) ...
                    + FD*(-2*FD*(cos(f).^2).*cosh(n).*sinh(n).*(kf.^2) + 1i*cos(f).*sinh(n).*kff + 2*1i*cos(f).*cosh(n).*kn ...
                    - 4*FD*(cos(f).^2).*cosh(n).*sinh(n).*(kn.^2) - 2*1i*cosh(n).*sin(f).*knf + 2*kf.*(-1i*sin(f).*sinh(n) ...
                    + FD*cos(f).*(cosh(n).^2).*(sin(f).*kn - cos(f).*knf)) + 1i*cos(f).*cosh(n).*knff + 3*1i*cos(f).*sinh(n).*knn ...
                    - 2*FD*((cos(f).*cosh(n)).^2).*kn.*knn + 1i*cos(f).*cosh(n).*k3n) - 2*h.*hn.*s)./(h.^2) ;
            end
            if nargout > 2
                snn = (-2*(FD^2)*(cos(f).^2).*cosh(2*n).*kf.^2 + 1i*FD*cos(f).*cosh(n).*kff + 2*1i*FD*cos(f).*sinh(n).*kn ...
                    - (6*(FD*cos(f).*cosh(n)).^2 - 2*h.^2 + 2*(FD*cosh(n).*sin(f)).^2 + 8*(FD*cos(f).*sinh(n)).^2).*kn.^2 ...
                    - 4*1i*FD*sin(f).*sinh(n).*knf + 4*(FD^2)*cos(f).*(cosh(n).^2).*sin(f).*kn.*knf - 2*(FD*cos(f).*cosh(n).*knf).^2 ...
                    + 2*1i*FD*cos(f).*sinh(n).*knff + 2*(k.^2).*(-(FD^2)*cosh(2*n) + hn.^2 + h.*hnn) + 5*1i*FD*cos(f).*cosh(n).*knn ...
                    - 7*((FD*cos(f)).^2).*sinh(2*n).*kn.*knn - 2*(FD*cos(f).*cosh(n).*knn).^2 - 2*1i*FD*cosh(n).*sin(f).*knnf ...
                    + FD*cosh(n).*kf.*(-2*1i*sin(f) + 2*FD*cos(f).*(4*sin(f).*sinh(n).*kn + cosh(n).*sin(f).*knn - cos(f).*(4*sinh(n).*knf + cosh(n).*knnf))) ...
                    + 1i*FD*cos(f).*cosh(n).*knnff + 4*1i*FD*cos(f).*sinh(n).*k3n - 2*((FD*cos(f).*cosh(n)).^2).*kn.*k3n ...
                    + 0.5*k.*(4*(FD^2).*cosh(2*n).*sin(2*f).*kf - 4*((FD^2)*(3 + cos(2*f)).*sinh(2*n) - 4*h.*hn).*kn - 6*(FD^2)*cosh(2*n).*knn ...
                    + 4*(h.^2).*knn + 2*(FD^2)*cosh(n).*sin(2*f).*(4*sinh(n).*knf + cosh(n).*knnf) - (FD^2)*sinh(2*n).*k3n ...
                    - (FD^2)*cos(2*f).*((-2 + 4*cosh(2*n)).*knn + sinh(2*n).* k3n)) + 1i*FD*cos(f).*cosh(n).*k4n ...
                    - (4*h.*hn.*sn + 2*h.*hnn.*s + 2*(hn.^2).*s) )./(h.^2);
            end
            if nargout > 3
                sf  = ((k.^2).*(-2*(FD^2)*cos(f).*sin(f) + 2*h.*hf) + k.*(2*(h.^2 + (FD^2)*(cos(2*f) - cosh(2*n).*sin(f).^2)).*kf ...
                    + (FD^2)*(cosh(n).*sin(2*f).*(cosh(n).*kff + 2*sinh(n).*kn) - (cos(f).^2).*sinh(2*n).*knf)) + FD*(4*FD*cos(f).*(cosh(n).^2).*sin(f).*kf.^2 ...
                    - 3*1i*cosh(n).*sin(f).*kff + 1i*cos(f).*cosh(n).*k3f - 2*1i*sin(f).*sinh(n).*kn + FD*(cosh(n).^2).*sin(2*f).*kn.^2 ...
                    - 2*cos(f).*cosh(n).*kf.*(1i + FD*cos(f).*(cosh(n).*kff + sinh(n).*kn)) + 2*1i*cos(f).*sinh(n).*knf - 2*FD*((cos(f).*cosh(n)).^2).*kn.*knf ...
                    - 1i*cosh(n).*sin(f).*knn + 1i*cos(f).*cosh(n).*knnf) - 2*h.*hf.*s)./(h.^2);
            end
            if nargout > 4
                sff = (((FD^2)*(4*cos(2*f) + (-1 + 3*cos(2*f)).*cosh(2*n)) + 2*h.^2).*kf.^2 + 2*(k.^2).*(-(FD^2)*cos(2*f) + hf.^2 + h.*hff) ...
                    + FD*cosh(n).*kf.*(2*1i*sin(f) + 7*FD*cosh(n).*sin(2*f).*kff - 2*FD*cos(f).*(cos(f).*cosh(n).*k3f ...
                    + 2*sinh(n).*(-2*sin(f).*kn + cos(f).*knf))) + 0.5*k.*(4*(-(FD^2)*(3 + cosh(2*n)).*sin(2*f) + 4*h.*hf).* kf ...
                    + (2*(FD^2)*(3*cos(2*f) + (-1 + 2*cos(2*f)).*cosh(2*n)) + 4*h.^2).*kff + 2*(FD^2)*(cosh(n).*sin(2*f).*(cosh(n).*k3f + 4*sinh(n).* knf) ...
                    + sinh(2*n).*(2*cos(2*f).*kn - (cos(f).^2).*knff))) + FD*(-2*FD*(cos(f).*cosh(n).*kff).^2 - cos(f).*cosh(n).*kff.*(5*1i + 2*FD*cos(f).*sinh(n).*kn) ...
                    + 4*FD*(cosh(n).^2).*sin(2*f).*kn.*knf + FD*cos(2*f).*(cosh(n).^2).*(2*kn.^2 - knf.^2 - kn.*knff) - FD*(cosh(n).^2).*(knf.^2 + kn.*knff) ...
                    - 2*1i*sin(f).*(2*sinh(n).*knf + cosh(n).*(2*k3f + knnf)) + 1i*cos(f).*(2*sinh(n).*(-kn + knff) + cosh(n).*(k4f - knn + knnff))) ...
                    - (4*h.*hf.*sf + 2*h.*hff.*s + 2*(hf.^2).*s) )./(h.^2);
            end
            
        end
    end
    
end
