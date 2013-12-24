classdef StarShapedMetrics < AbstractMetrics
		% unfinished for future use

    properties
        TheBody;
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
        
        function obj = StarShapedMetrics(Body,theta)
            obj.TheBody = Body;
            
             [curv,ds_curv,dsds_curv] = obj.Curvature(theta);
            
%             obj.h   = 
%             obj.hn  = 
%             obj.hnn = 
%             obj.h3n = 
%             obj.h4n = 
%                     
%             
%             obj.hf  = 
%             obj.hff = 
%             
%             obj.h3f = 
%             obj.h4f = 
%                     
            
            
		end        
         end
end