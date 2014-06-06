classdef Parametrization < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Public)
        X_FuncWithDeriatives;
        Y_FuncWithDeriatives;
    end
    
    methods
        function obj = ParametricEllipse(XHandle, YHandle)
            obj.X_FuncWithDeriatives = XHandle;
            obj.Y_FuncWithDeriatives = YHandle;
        end
    end
    
end

