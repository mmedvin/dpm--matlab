classdef ConstantWaveNumber < Tools.Common.FunctionWithDerivatives
    properties
        k;        
    end
    
    methods(Static = true)
        function [k,kr] = kkr(r,r0,k0)
            k=obj.k;
            kr=0;
        end
    end
        
    
    methods
        
        function [k,kr,krr,k3r,k4r,k5r] = Derivatives(obj)
            k   = obj.k;
            kr  = 0;
            krr = 0;
            k3r = 0;
            k4r = 0;
            k5r = 0;
        end
        
        function obj=ConstantWaveNumber(~,k) %(k0,r,r0)
            obj.k=k;
        end
        
    end
    
    
end
