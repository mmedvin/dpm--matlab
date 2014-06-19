classdef ParametricHeart < Tools.Parameterizations.AbstractParameterization
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        XHandle;
        YHandle;
    end
    
    methods
        function obj = ParametricHeart(Params)
            obj.XHandle = Tools.Parameterizations.HeartX(Params.e,Params.p);
            obj.YHandle = Tools.Parameterizations.HeartY(Params);
        end
    end
end
