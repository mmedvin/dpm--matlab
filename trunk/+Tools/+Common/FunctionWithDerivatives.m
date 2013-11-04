classdef FunctionWithDerivatives < handle    
    methods (Abstract)
        Derivatives(obj);
    end
end