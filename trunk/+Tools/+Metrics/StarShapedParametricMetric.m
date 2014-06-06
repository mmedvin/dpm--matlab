classdef StarShapedParametricMetric < Tools.Metrics.AbstractMetrics
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = protected)
        Xt;
        Yt;
    end
    
    methods
        function obj = StarShapedParametricMetric(XtFuncHndl,YtFuncHndl)
            obj.Xt = XtFuncHndl;
            obj.Yt = YtFuncHndl;
        end
        
        function [h,ht,htt,h3t] = metrics(obj,t)
            [~,xt,xtt,x3t,x4t] = obj.Xt.Derivatives(t);
            [~,yt,ytt,y3t,y4t] = obj.Yt.Derivatives(t);
            
            h   = sqrt(xt.^2 + yt.^2);
            ht  = (xt.*xtt + yt.*ytt)./h;
            htt = (xtt.^2 + ytt.^2 + xt.*x3t + yt.*y3t - ht.^2)./h;
            h3t = (3*xtt.*x3t + 3*ytt.*y3t + xt.*x4t + yt.*y4t +  - 3*ht.*htt)./h;
        end
        
        function [curv,ds_curv,dsds_curv] = curvature(obj,t)
            [~,xt,xtt,x3t,x4t] = obj.Xt.Derivatives(t);
            [~,yt,ytt,y3t,y4t] = obj.Yt.Derivatives(t);
            [h,ht,htt]     = obj.metrics(t);
            
            h3 = (h.^3);
            curv = (xt.*ytt - yt.*xtt)./h3;
            ct   = -3*curv.*ht./h + (xt.*y3t - yt.*x3t)./h3;
            ctt  = -3*curv.*(2*(ht.^2)./h + htt)./h - 6*ht .*ct./h + (xtt.*y3t + xt.*y4t - x3t.*ytt - yt.*x4t)./h3;
            
            ds_curv     = ct./h;
            dsds_curv   = (ctt - ht.*ds_curv)./(h.^2);     
            
            curv  = -curv;
            ds_curv = -ds_curv;
            dsds_curv = -dsds_curv;
            
        end
        
        
    end
    
end

