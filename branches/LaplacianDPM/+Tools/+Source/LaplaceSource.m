classdef LaplaceSource < Tools.Source.SuperSource
    properties (Dependent = true)
        Source;%Fn;Ff;Fnn;Fff;    % WN;
    end
    
    properties(Access = protected)       
        Scatterer; 
    end
    
    methods
        function [F,Fn,Ff,Fnn,Fff] = Derivatives(obj)
            
            if obj.IsDummy
                F=0;
                Fn=0;
                Ff=0;
                Fnn=0;
                Fff=0;
            else
                r = scatterer.r;
                F   = 8*r.^2 + 4;
                
                if nargout>1, Fn  = 16*r;  end
                if nargout>2, Ff  = 0   ;  end
                if nargout>3, Fnn = 16  ;  end
                if nargout>4, Fff = 0   ;  end
            end
            
        end
            
        function obj = LaplaceSource(Scatterer,~)%,Coeffs)
            obj = obj@Tools.Source.SuperSource(Scatterer,[]);                      
        end
        
        
    end
    
end

