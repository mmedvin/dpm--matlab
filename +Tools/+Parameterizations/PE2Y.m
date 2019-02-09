classdef PE2Y < Tools.Common.FunctionWithDerivatives
%XHandle for Parametric Ellipse 2    
    properties
        C;%acost
        S;%bsint
        ycenter;
        rot;
        a;
        b;
    end
    
        methods        
        function obj = PE2Y(a,b,ycenter,rotation)
            
            obj.C=Tools.Parameterizations.AcosBtWD(a,1);
            obj.S=Tools.Parameterizations.AsinBtWD(b,1);
            obj.ycenter = ycenter;
            obj.rot = rotation;
            obj.a=a;
            obj.b=b;
            
        end
        
        function varargout = Derivatives(obj,t)
            nout = nargout;
            [Cs{1:nout }]= obj.C.Derivatives(t);
            [Sn{1:nout }] = obj.S.Derivatives(t);
            
            varargout(1)= {obj.ycenter + Cs{1}.*sin(obj.rot) + Sn{1}.*cos(obj.rot)};
            
            for n=2:nout 
                 varargout(n)= {Cs{n}.*sin(obj.rot) + Sn{n}.*cos(obj.rot)};
            end
            
            %            nout = nargout;
            
%             if nout > 0, varargout(1)={           obj.a*cos(obj.b*t) + obj.c*sin(obj.d*t) };  end %u
%             if nout > 1, varargout(2)={-(obj.b  )*obj.a*sin(obj.b*t) +(obj.d  )*obj.c*cos(obj.d*t) };  end %ut
%             if nout > 2, varargout(3)={-(obj.b^2)*varargout{1}        };  end %utt
%             if nout > 3, varargout(4)={-(obj.b^2)*varargout{2}        };  end %u3t
%             if nout > 4, varargout(5)={ (obj.b^4)*varargout{1}        };  end %u4t
            
            % obj.u5t     = b^4*obj.ut;
            % obj.u6t     = b^4*obj.utt;
            % obj.u7t     = b^4*obj.u3t;
        end
        
    end

    
    
end