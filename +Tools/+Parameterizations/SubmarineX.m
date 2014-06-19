classdef SubmarineX < Tools.Common.FunctionWithDerivatives
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        a;c;p;d;e;
        TowerWD
    end
    
    methods
        function obj = SubmarineX(Params)
            obj.a=Params.a;
            obj.c=Params.c;
            obj.p=Params.p;
            obj.d=Params.d;
            obj.e=Params.e;
            
            obj.TowerWD = Tools.Common.FunctionWithDerivatives.SubmarineTowerPatch(Params);
        end
        function varargout = Derivatives(obj,t)
            nout = nargout;
            
           TD{1:nout+1} = obj.TowerWD.Derivatives(t);%TD{n} is n-1 derivative of TD 
            
            if nout > 0, varargout(1)={ obj.a*cos(t).*TD{1}                                                                         };  end %u
            if nout > 1, varargout(2)={-obj.a*sin(t).*TD{1} +   obj.a*cos(t).*TD{2}                                                 };  end %ut
            if nout > 2, varargout(3)={-obj.a*cos(t).*TD{1} - 2*obj.a*sin(t).*TD{2} +   obj.a*cos(t).*TD{3}                        	};  end %utt
            if nout > 3, varargout(4)={ obj.a*sin(t).*TD{1} - 3*obj.a*cos(t).*TD{2} - 3*obj.a*sin(t).*TD{3} + obj.a*cos(t).*TD{4}   };  end %u3t
            if nout > 4, varargout(5)={ obj.a*cos(t).*TD{1} + 4*obj.a*sin(t).*TD{2} - 6*obj.a*cos(t).*TD{3} ...
                                                            - 4*obj.a*sin(t).*TD{4} +   obj.a*cos(t).*TD{5}                     	};  end %u4t
            
            % obj.u5t     = b^4*obj.ut;
            % obj.u6t     = b^4*obj.utt;
            % obj.u7t     = b^4*obj.u3t;
        end
    end
    
end

