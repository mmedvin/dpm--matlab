classdef NestedPolarScatterer < Tools.Scatterer.NestedScatterer
    
    properties
        r0;
        r1;

        
        r;
        th;
        R;
        Th;
    end
    
    methods
        function obj =  NestedPolarScatterer(Grid,Params)
            
            Arguments.InteriorScatterer  = @Tools.Scatterer.PolarScatterer;
            Arguments.ExteriorScatterer  = @Tools.Scatterer.PolarScatterer;
            Arguments.IntScatParams      = struct('r0',Params.r0, 'Stencil',Params.Stencil);
            Arguments.ExtScatParams      = struct('r0',Params.r1, 'Stencil',Params.Stencil);
            Arguments.Grid               = Grid;
                        
            obj      = obj@Tools.Scatterer.NestedScatterer(Arguments);            
            obj.r0   = Params.r0;
            obj.r1   = Params.r1;
        end
%         function L = get.Inside(obj)
%             L = (obj.R >= obj.r0) && (obj.R <= obj.r1);
%         end
%         function L = get.Outside(obj)
%             L = (obj.R < obj.r0) && (obj.R > obj.r1);
%         end
        
        function Expansion(obj,Xi0,Xi1,F,WaveNumber)
            keyboard; %error('tbd');%
        end
        
        function val = get.r(obj)
            val = {obj.InteriorScatterer.r, obj.ExteriorScatterer.r};
        end
        
        function val = get.th(obj)
            val = {obj.InteriorScatterer.th, obj.ExteriorScatterer.th};
        end
        
        function val = get.R(obj)
            val = obj.Grid.R;
        end
        
        
        function val = get.Th(obj)
            val = obj.Grid.Theta();
        end
    end
end