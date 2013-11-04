classdef SuperExact < Tools.Common.FunctionWithDerivatives
    properties
        Scatterer; 
        WaveNumber;
    end
    methods
        function obj = SuperExact(Scatterer, WaveNumber)
            obj.Scatterer  = Scatterer;
            obj.WaveNumber = WaveNumber;
        end
    end
end
