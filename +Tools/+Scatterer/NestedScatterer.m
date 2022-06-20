classdef NestedScatterer < Tools.Scatterer.SuperScatterer
    
     properties(SetAccess=protected)
          Grid;

        
        Mp;
        Mm;
        Np;
        Nm;
        
        Inside;
        Outside;
        
        GridGamma;
                
        In; %can't find better name yet...
        Out;
        
        %%%%%%%%%%%%%
        BasisArg;
        %TheScatterer;  
        MetricsAtScatterer=0;
        
        InteriorScatterer;
        ExteriorScatterer;
    end
     
%     methods(Access=protected)
%     
%     end
    
    
    methods
        function obj = NestedScatterer(Arguments)            
            obj.Grid = Arguments.Grid;
            obj.InteriorScatterer = Arguments.InteriorScatterer(Arguments.Grid,Arguments.IntScatParams);
            obj.ExteriorScatterer = Arguments.ExteriorScatterer(Arguments.Grid,Arguments.ExtScatParams);

            %obj.InteriorScatterer.r=-obj.InteriorScatterer.r;
            %obj.ExteriorScatterer.r=-obj.ExteriorScatterer.r;
            
            obj.BasisArg        = {obj.InteriorScatterer.BasisArg    , obj.ExteriorScatterer.BasisArg       };
            % obj.TheScatterer    = {obj.InteriorScatterer.TheScatterer, obj.ExteriorScatterer.TheScatterer   };
            obj.GridGamma       = union(obj.InteriorScatterer.GridGamma   , obj.ExteriorScatterer.GridGamma );
            
            obj.Inside  = xor(obj.InteriorScatterer.Inside , obj.ExteriorScatterer.Inside  );
            obj.Outside = ~obj.Inside;
            
            obj.Mp = find(obj.Inside());
            obj.Mm = find(obj.Outside());
            obj.Np = intersect(obj.InteriorScatterer.Nm, obj.ExteriorScatterer.Np);
            obj.Nm = union    (obj.InteriorScatterer.Np, obj.ExteriorScatterer.Nm);

        end
        
        function in = get.In(obj)            
            in = setdiff(setdiff(obj.Mp,obj.GridGamma{1}),obj.GridGamma{1});            
        end
        
        function out = get.Out(obj)
            out = setdiff(setdiff(obj.Mm,obj.GridGamma{1}),obj.GridGamma{1});
        end                
    end    
end
