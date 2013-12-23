classdef SuperHelmholtzSource < Tools.Common.FunctionWithDerivatives
       
    properties(Access = protected)  
        IsDummy;
    end
    
    methods
        function obj = SuperHelmholtzSource()
            obj.IsDummy = true;
        end
        function varargout = Derivatives(obj)
            assert(obj.IsDummy);
            
            nout = nargout;
            varargout(1:nout)={0};    
        end
    end
end

