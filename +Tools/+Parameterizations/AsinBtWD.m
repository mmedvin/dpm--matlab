classdef AsinBtWD < Tools.Common.FunctionWithDerivatives
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %coefficients of a*sin(b*t)
        a;
        b;
    end
    
    methods        
        function obj = AsinBtWD(a,b)
            obj.a=a;
            obj.b=b;

        end
        
        function varargout = Derivatives(obj,t)
            nout = nargout;
            
            if nout > 0, varargout(1)={           obj.a*sin(obj.b*t)  };  end %u
            if nout > 1, varargout(2)={ (obj.b  )*obj.a*cos(obj.b*t)  };  end %ut
            if nout > 2, varargout(3)={-(obj.b^2)*varargout{1}        };  end %utt
            if nout > 3, varargout(4)={-(obj.b^2)*varargout{2}        };  end %u3t
            if nout > 4, varargout(5)={ (obj.b^4)*varargout{1}        };  end %u4t
            % obj.u5t     = b^4*obj.ut;
            % obj.u6t     = b^4*obj.utt;
            % obj.u7t     = b^4*obj.u3t;
        end
        
    end
    
end
