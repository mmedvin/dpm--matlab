classdef ParametricStar < Tools.Parameterizations.AbstractParameterization

    
    properties(Access = public)
        XHandle;
        YHandle;
    end
    
    methods
        function obj = ParametricStar()
            obj.XHandle = Tools.Parameterizations.StarX();
            
            obj.YHandle = Tools.Parameterizations.StarY();
        end
    end
end