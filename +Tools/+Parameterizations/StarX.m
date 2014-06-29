classdef StarX < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function varargout = Derivatives(obj,t)
            nout = nargout;
            
            a=7/6;
            b=1/6;
            c=4;
            if nout > 0, varargout(1)={	a*cos(t) + b*cos(c*t)           };  end %u
            if nout > 1, varargout(2)={-a*sin(t) - c*b*sin(c*t)         };  end %ut
            if nout > 2, varargout(3)={-a*cos(t) - (c^2)*b*cos(c*t)     };  end %utt
            if nout > 3, varargout(4)={ a*sin(t) + (c^3)*b*sin(c*t)     };  end %u3t
            if nout > 4, varargout(5)={ a*cos(t) + (c^4)*b*cos(c*t)     };  end %u4t
            
            % obj.u5t     = b^4*obj.ut;
            % obj.u6t     = b^4*obj.utt;
            % obj.u7t     = b^4*obj.u3t;
        end
    end
    
end

%1/6 {5 Cos[t] + Cos[4 t], 5 Sin[t] - Sin[4 t]}