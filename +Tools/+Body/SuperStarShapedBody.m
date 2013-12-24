classdef SuperStarShapedBody < Tools.Common.FunctionWithDerivatives 
    properties(Access = protected)
        r;
        dr;
        drr;
        d3r;
        d4r;
        d5r;
    end
    methods
        function [R,DR,DRR,D3R,D4R,D5R] = Derivatives(obj)
            R   = obj.r;
            DR  = obj.dr;
            DRR = obj.drr;
            D3R = obj.d3r;
            D4R = obj.d4r;
            D5R = obj.d5r;
        end
    end
end
    
