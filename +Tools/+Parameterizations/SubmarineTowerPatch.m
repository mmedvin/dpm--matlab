classdef SubmarineTowerPatch < Tools.Common.FunctionWithDerivatives
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c;
        p;
        %d;e; 
    end
    
    methods
        
        function obj = SubmarineTowerPatch(c,p)
            obj.c=c;
            obj.p=p;           
        end
        
        function varargout = Derivatives(obj,t)
            nout = nargout;
            
            if nout > 0, varargout(1)={  1 + obj.c.*( (sin(t/2) + cos(t/2))/sqrt(2)).^obj.p                                                  };  end %u
            
            if nout > 1, varargout(2)={  (2^(1/2).*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 1))/2 	};  end %ut
       
            if nout > 2
                varargout(3)={  (obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).^2.*(obj.p - 1).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 2))/2 ...
                              - (2^(1/2).*obj.c.*obj.p.*(cos(t/2)/4 + sin(t/2)/4).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 1))/2                                 
                };  
            end %utt
            if nout > 3
                varargout(4)={  (2^(1/2).*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).^3.*(obj.p - 1).*(obj.p - 2).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 3))/4 ...
                              - (3.*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).*(cos(t/2)/4 + sin(t/2)/4).*(obj.p - 1).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 2))/2 ...
                              - (2^(1/2).*obj.c.*obj.p.*(cos(t/2)/8 - sin(t/2)/8).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 1))/2 
                    };  
            end %u3t
            
            if nout > 4
                arr  = (3.*obj.c.*obj.p.*(cos(t/2)/4 + sin(t/2)/4).^2.*(obj.p - 1).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 2))/2 ...
                              + (2^(1/2).*obj.c.*obj.p.*(cos(t/2)/16 + sin(t/2)/16).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 1))/2 ...
                              - 2.*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).*(cos(t/2)/8 - sin(t/2)/8).*(obj.p - 1).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 2) ...
                              +(obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).^4.*(obj.p - 1).*(obj.p - 2).*(obj.p - 3).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 4))/4 ...
                              - (3.*2^(1/2).*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).^2.*(cos(t/2)/4 ...
                              + sin(t/2)/4).*(obj.p - 1).*(obj.p - 2).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 3))/2;
                          
                varargout(5)={arr                    };
            end %u4t
            if nout > 5
                arr =  (2^(1/2).*obj.c.*obj.p.*(cos(t/2)/32 - sin(t/2)/32).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 1))/2 ...
                              + 5.*obj.c.*obj.p.*(cos(t/2)/4 + sin(t/2)/4).*(cos(t/2)/8 - sin(t/2)/8).*(obj.p - 1).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 2) ... 
                              +(5.*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).*(cos(t/2)/16 + sin(t/2)/16).*(obj.p - 1).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 2))/2 ...
                              -(5.*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).^3.*(cos(t/2)/4 + sin(t/2)/4).*(obj.p - 1).*(obj.p - 2).*(obj.p - 3).*((2^(1/2).*(cos(t/2) ...
                              + sin(t/2)))/2).^(obj.p - 4))/2 + (15.*2^(1/2).*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).*(cos(t/2)/4 ... 
                              + sin(t/2)/4).^2.*(obj.p - 1).*(obj.p - 2).*((2^(1/2).*(cos(t/2) + sin(t/2)))/2).^(obj.p - 3))/4 ...
                              -(5.*2^(1/2).*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).^2.*(cos(t/2)/8 - sin(t/2)/8).*(obj.p - 1).*(obj.p - 2).*((2^(1/2).*(cos(t/2) ...
                              + sin(t/2)))/2).^(obj.p - 3))/2 + (2^(1/2).*obj.c.*obj.p.*(cos(t/2)/2 - sin(t/2)/2).^5.*(obj.p - 1).*(obj.p - 2).*(obj.p - 3).*(obj.p - 4).*((2^(1/2).*(cos(t/2) ...
                              + sin(t/2)))/2).^(obj.p - 5))/8 ;
                varargout(6)={ arr };  
            end %u5t
            
            % obj.u5t     = b^4.*obj.ut;
            % obj.u6t     = b^4.*obj.utt;
            % obj.u7t     = b^4.*obj.u3t;
        end
    end
    
end




%  methods
%         function varargout = Derivatives(obj,t)
%             nout = nargout;
%             
%             StepF  = 
%             
%             
%             if nout > 0, varargout(1)={  (  1 + (obj.c.* (sin(t).^obj.p) )./(obj.d .* exp(obj.e .* (t - pi)) + 1)  ) };  end %u
%             
%             if nout > 1, varargout(2)={ ...
%                 (obj.c.*obj.p.* cos(t)..* (sin(t).^(obj.p-1)) )./(1 + obj.d.*exp(obj.e.*(t-pi))) ...
%                -(obj.c.*obj.d.*obj.e.*exp(obj.e.*(t-pi))..* (sin(t).^obj.p) )./(1 + obj.d.*exp(obj.e (t-pi))).^2;           };  end %ut
%        
%             if nout > 2, varargout(3)={ ...      
%                     -((2.*obj.c.*obj.d.*obj.e.*exp(obj.e.* (t-pi))..* obj.p.* cos(t)..*(sin(t).^(obj.p-1)))./((1 + obj.d.*exp(obj.e.* (t-pi))).^2)) ...
%                     + obj.c.*((2.* ((obj.d.*obj.e)^2).*exp(2.*obj.e.*(t-pi)))./((1 + obj.d.*exp(obj.e.* (t-pi)))^3) ...
%                     - (obj.d.* (obj.e^2).*exp(obj.e.* (t-pi)))./(1 + obj.d.*exp(obj.e.* (t-pi)))^2) Sin[t]^
%   p + (c ((-1 + p) p Cos[t]^2 Sin[t]^(-2 + p) - p Sin[t]^p))/(
%  1 + d E^(e (-\[Pi] + t)))
%                 
%                 
%                 };  end %utt
%             if nout > 3, varargout(4)={       };  end %u3t
%             if nout > 4, varargout(5)={        };  end %u4t
%             if nout > 5, varargout(6)={        };  end %u5t
%             
%             % obj.u5t     = b^4.*obj.ut;
%             % obj.u6t     = b^4.*obj.utt;
%             % obj.u7t     = b^4.*obj.u3t;
%         end
%     end


