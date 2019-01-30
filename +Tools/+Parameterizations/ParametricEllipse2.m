classdef ParametricEllipse2 < Tools.Parameterizations.AbstractParameterization
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        XHandle;
        YHandle;
    end
    
    methods
        function obj = ParametricEllipse2(Params)
            obj.XHandle = Tools.Parameterizations.PE2X(Params.a,Params.b,Params.xcenter,Params.rotation);
            obj.YHandle = Tools.Parameterizations.PE2Y(Params.a,Params.b,Params.ycenter,Params.rotation);
        end
    end
    
end
