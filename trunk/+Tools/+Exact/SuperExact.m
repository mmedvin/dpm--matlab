classdef SuperExact < Tools.Common.FunctionWithDerivatives
    properties
        Scatterer; 
        Coeffs;
    end
    methods
        function obj = SuperExact(Scatterer, Coeffs)
            obj.Scatterer  = Scatterer;
            obj.Coeffs = Coeffs;
        end
    end
end
