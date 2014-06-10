classdef ParametricKite < Tools.Parameterizations.AbstractParameterization
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        XHandle;
        YHandle;
    end
    
    methods
        function obj = ParametricKite(Params)
            obj.XHandle = Tools.Parameterizations.AcosPlusHalfBCos2tMinusHalfBWD(Params.a,Params.b);
            obj.YHandle = Tools.Parameterizations.AsinBtWD(Params.c,1);
        end
    end
end
