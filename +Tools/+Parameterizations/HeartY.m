classdef HeartY < Tools.Common.FunctionWithDerivatives
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %coefficients of a*cos(b*t)
        a;
        b;
        c;
        d;        
    end
    
    methods        
        function obj = HeartY(Params)
            obj.a=Params.a;
            obj.b=Params.b;
            obj.c=Params.c;
            obj.d=Params.d;
        end
        
        function varargout = Derivatives(obj,t)
            nout = nargout;
            
            %13 Cos[t] - 5 Cos[2 t] - 2 Cos[3 t] - Cos[4 t]
            
            if nout > 0, varargout(1)={    obj.a*cos(t) + obj.b*cos(2*t) +  obj.c*cos(3*t) + obj.d*cos(4*t)                 };  end %u
            if nout > 1, varargout(2)={- 2*obj.b*sin(2*t) - 3*obj.c*sin(3*t) - 4*obj.d*sin(4*t) - obj.a*sin(t)              };  end %ut
            if nout > 2, varargout(3)={- obj.a*cos(t) - 4*obj.b*cos(2*t) - 9*obj.c*cos(3*t) - 16*obj.d*cos(4*t)             };  end %utt
            if nout > 3, varargout(4)={  8*obj.b*sin(2*t) + 27*obj.c*sin(3*t) + 64*obj.d*sin(4*t) + obj.a*sin(t)            };  end %u3t
            if nout > 4, varargout(5)={  obj.a*cos(t) + 16*obj.b*cos(2*t) + 81*obj.c*cos(3*t) + 256*obj.d*cos(4*t)          };  end %u4t
            
            % obj.u5t     = b^4*obj.ut;
            % obj.u6t     = b^4*obj.utt;
            % obj.u7t     = b^4*obj.u3t;
        end
        
    end
    
end

