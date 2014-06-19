classdef HeartX < Tools.Common.FunctionWithDerivatives
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %coefficients of a*sin(b*t)
        e;
        p;
    end
    
    methods        
        function obj = HeartX(e,p)
            obj.e=e;
            obj.p=p;

        end
        
        function varargout = Derivatives(obj,t)
            nout = nargout;
            
            if nout > 0, varargout(1)={  obj.e*sin(t).^obj.p                    };  end %u
            if nout > 1, varargout(2)={  obj.e*obj.p*cos(t).*sin(t).^(obj.p - 1)    };  end %ut
            if nout > 2, varargout(3)={  obj.e*obj.p*cos(t).^2.*sin(t).^(obj.p - 2).*(obj.p - 1) - obj.e*obj.p.*sin(t).*sin(t).^(obj.p - 1)      };  end %utt
            if nout > 3, 
                arr = obj.e*obj.p*cos(t).^3.*sin(t).^(obj.p - 3).*(obj.p - 1).*(obj.p - 2) - obj.e*obj.p*cos(t).*sin(t).^(obj.p - 1) ...
                                        -3*obj.e*obj.p*cos(t).*sin(t).*sin(t).^(obj.p - 2).*(obj.p - 1);
                varargout(4)={  arr      };  
            end %u3t
            if nout > 4, 
                arr = obj.e*obj.p*sin(t).*sin(t).^(obj.p - 1) - 4*obj.e*obj.p*cos(t).^2.*sin(t).^(obj.p - 2).*(obj.p - 1) ...
                                       + 3*obj.e*obj.p*sin(t).^(obj.p - 2).*sin(t).^2*(obj.p - 1) ...
                                       - 6*obj.e*obj.p*cos(t).^2.*sin(t).*sin(t).^(obj.p - 3)*(obj.p - 1).*(obj.p - 2) ...
                                       + obj.e*obj.p*cos(t).^4.*sin(t).^(obj.p - 4)*(obj.p - 1)*(obj.p - 2)*(obj.p - 3);
                varargout(5)={  arr       };  
            end %u4t
            % obj.u5t     = b^4*obj.ut;
            % obj.u6t     = b^4*obj.utt;
            % obj.u7t     = b^4*obj.u3t;
        end
        
    end
    
end
