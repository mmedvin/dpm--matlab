classdef NestedExtension < Tools.Extensions.SuperExtension 
 
    properties
        InteriorExtension;
        ExteriorExtension;
    end
    methods
        function obj = NestedExtension(Arguments)
            
%             obj.Basis     = Arguments.Basis;
%             obj.Coeffs    = Arguments.Coeffs;
%             obj.Scatterer = Arguments.Scatterer;
%             obj.Grid      = Arguments.Grid;
            
            IntArguments            = Arguments;
            IntArguments.Scatterer  = Arguments.Scatterer.InteriorScatterer;
            IntArguments.Basis      = Arguments.Basis.Interior;
            
            obj.InteriorExtension   = Arguments.IntExtension(IntArguments);
            
            ExtArguments            = Arguments;
            ExtArguments.Scatterer  = Arguments.Scatterer.ExteriorScatterer;
            ExtArguments.Basis      = Arguments.Basis.Exterior;

            obj.ExteriorExtension   = Arguments.ExtExtension(ExtArguments);

        end  
        
        function ExpandSource(obj,SourceHandle,SourceParams)
            %keyboard;
            obj.InteriorExtension.ExpandSource(SourceHandle,SourceParams);
            obj.ExteriorExtension.ExpandSource(SourceHandle,SourceParams);
            obj.Wf= [  obj.InteriorExtension.Wf; obj.ExteriorExtension.Wf ];
        end
        
        function Expand(obj)
            obj.InteriorExtension.Expand();
            obj.ExteriorExtension.Expand();
            
            obj.W = {obj.InteriorExtension.W{1}, ...
                     obj.InteriorExtension.W{2}; ...
                     obj.ExteriorExtension.W{1}, ...
                     obj.ExteriorExtension.W{2}
                    };           
        end
        
        function val = Expansion(obj,xi0,xi1,Src)
            keyboard; 
        end
        
        function [xi0j,xi1j] = ExpandedBasis(obj,n)
            keyboard;
        end
    end
    
end

