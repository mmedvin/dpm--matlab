classdef NavierStokesCombinedExtension < Tools.Extensions.SuperExtension 
    %NSPSIEXTENSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        OmegaExtension;
        PsiExtension;
    end
    
    methods
        function obj = NavierStokesCombinedExtension(Arguments)
            
            OmegaArguments            = Arguments;
            OmegaArguments.Scatterer  = Arguments.Scatterer;
            OmegaArguments.Basis      = Arguments.Basis;
            
            obj.OmegaExtension   = Arguments.OmegaExtension(OmegaArguments);
            
            PsiArguments            = Arguments;
            PsiArguments.Scatterer  = Arguments.Scatterer;
            PsiArguments.Basis      = Arguments.Basis;
            
            obj.PsiExtension   = Arguments.PsiExtension(PsiArguments);
            
        end
        
        function ExpandSource(obj,SourceHandle,SourceParams)
            %keyboard;
            obj.OmegaExtension.ExpandSource(SourceHandle,SourceParams);
            obj.PsiExtension.ExpandSource(SourceHandle,SourceParams);
            obj.Wf= [  obj.OmegaExtension.Wf; obj.PsiExtension.Wf ];
        end
        
        function Expand(obj)%???
            obj.OmegaExtension.Expand();
            obj.PsiExtension.Expand();
            
            obj.W = {obj.OmegaExtension.W{1}, ...
                     obj.OmegaExtension.W{2}; ...
                     obj.PsiExtension.W{1}, ...
                     obj.PsiExtension.W{2}
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

