classdef ConstantWaveNumber < Tools.Coeffs.AbstractCoeffs
    properties
        k;        
    end
    
    methods(Static = true)
        function [k,kr] = kkr(~,~,k0)
            k=k0;
            kr=0;
        end
    end
        
    
    methods                
        function varargout = Derivatives(obj,~)
            varargout(1)={obj.k};
            nout = nargout;
            varargout(2:nout)={0};
        end
        
        function obj=ConstantWaveNumber(Grid,Params)
            obj.k=Params.k;
            
            if isa(Grid, 'Tools.Grid.Grids') || isa(Grid, 'Tools.Scatterer.SingleScatterer')
                obj.k=Params.k*spones(ones(size(Grid.R)));
            end
        end
        
    end
    
    
end
