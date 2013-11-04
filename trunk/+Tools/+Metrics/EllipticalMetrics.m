classdef EllipticalMetrics < Tools.Metrics.AbstractMetrics
    properties
        h   =0;
        hn  =0;
        hnn =0;
        h3n =0;
        h4n =0;
        hf  =0;
        hff =0;
        h3f =0;
        h4f =0;
    end
    
    methods       
        function [h,hn,hnn,h3n,h4n,hf,hff,h3f,h4f] = metrics(obj)
            h   = obj.h;
            hn  = obj.hn;
            hnn = obj.hnn;
            h3n = obj.h3n;
            h4n = obj.h4n;
            hf  = obj.hf;
            hff = obj.hff;
            h3f = obj.h3f;
            h4f = obj.h4f;
        end
        
        function obj = EllipticalMetrics(FocalDistance,eta,phi)
            f = FocalDistance;
            
            obj.h   =  f*sqrt(sinh(eta).^2 + sin(phi).^2);
            obj.hn  = (f^2).*sinh(2*eta)./(2*obj.h);
            obj.hnn = (f^2).*cosh(2*eta)./obj.h - (obj.hn.^2)./obj.h;
            obj.h3n = (4 - (3* obj.hnn)./ obj.h) .* obj.hn;
            obj.h4n = (-((3* obj.h3n)./ obj.h) + (  3*obj.hnn.*obj.hn)./ ((obj.h).^2)).* obj.hn ...
                    + (4 - (3* obj.hnn)./obj.h).* obj.hnn;
            
            obj.hf  = (f^2).*sin(2*phi)./(2*obj.h);
            obj.hff = (f^2).*cos(2*phi)./obj.h - (obj.hf.^2)./obj.h;
            
            obj.h3f = (-4 - (3* obj.hff)./ obj.h) .* obj.hf;
            obj.h4f = (-((3* obj.h3f)./ obj.h) + (  3*obj.hff.*obj.hf)./ ((obj.h).^2)).* obj.hf ...
                    + (4 - (3* obj.hff)./obj.h).* obj.hff;

        end
    end
    
end
