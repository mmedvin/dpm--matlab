classdef SuperPolarTwoTupleExtension < Tools.Extensions.TwoTupleExtension
 % Equaion Based Polar Helmholtz 5 Order Extension
    
    
    properties(Access=protected)
        r0;
        dr;
        th;
    end
    
    methods
        function obj = SuperPolarTwoTupleExtension(Arguments)
            obj = obj@Tools.Extensions.TwoTupleExtension(Arguments);
        end

        function val = get.r0(obj)
            val = obj.Scatterer.r0;
        end
        
        function val = get.dr(obj)
            val = obj.Scatterer.dr;
        end
        
        function val = get.th(obj)
            val = obj.Scatterer.th;
        end
    end
    
end

