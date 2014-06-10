classdef ParametricEllipse < Tools.Parameterizations.AbstractParameterization
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        XHandle;
        YHandle;
    end
    
    methods
        function obj = ParametricEllipse(Params)
            obj.XHandle = Tools.Parameterizations.AcosBtWD(Params.a,1);
            obj.YHandle = Tools.Parameterizations.AsinBtWD(Params.b,1);
        end
    end
    
end
