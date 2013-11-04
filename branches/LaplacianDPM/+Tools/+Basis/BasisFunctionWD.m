classdef BasisFunctionWD < Tools.Common.FunctionWithDerivatives
    properties
        xi0;
        xi0t;
        xi0tt;
        xi0ttt;
        xi0tttt;
        xi0tttttt;
    end
    
    methods
        
        function [xi0,xi0f,xi0ff,xi0ffff,xi0ffffff,xi0fff] = Derivatives(obj)
            xi0     = obj.xi0;
            xi0f    = obj.xi0t;
            xi0ff   = obj.xi0tt;
            xi0ffff = obj.xi0tttt;
            xi0ffffff = obj.xi0tttttt;
            
            xi0fff   = obj.xi0ttt;
        end
        
    end
end
