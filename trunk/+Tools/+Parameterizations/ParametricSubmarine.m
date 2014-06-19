classdef ParametricSubmarine < Tools.Parameterizations.AbstractParameterization
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        XHandle;
        YHandle;
    end
    
    methods
        function obj = ParametricSubmarine(Params)
            obj.XHandle = Tools.Parameterizations.AcosBtWD(Params.a,1);
            
            obj.YHandle = Tools.Parameterizations.SubmarineY(Params);
        end
    end
end