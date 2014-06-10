classdef AcosPlusHalfBCos2tMinusHalfBWD < Tools.Common.FunctionWithDerivatives
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
   
    properties
        %coefficients of a*cos(b*t)
        a;
        b;
    end
   
%     properties(Access = protected)
%         ACost;
%         HalfBCos2t;
%        
%     end
   
    methods       
        function obj = AcosPlusHalfBCos2tMinusHalfBWD(a,b)
           
            %obj.ACost       = Tools.Parameterizations.AcosBtWD(a,1);
            %obj.HalfBCos2t  = Tools.Parameterizations.AcosBtWD(b/2,2);
           
            %a Cos[t] + 0.5 b (Cos[2 t] - 1)
           
            obj.a=a;
            obj.b=b;
        end
       
        function varargout = Derivatives(obj,t)
            nout = nargout;

            if nout > 0, varargout(1)={ obj.a*cos(t) + 0.5*obj.b*cos(2*t) - .5*obj.b    };  end %u
            if nout > 1, varargout(2)={-obj.a*sin(t) -     obj.b*sin(2*t)               };  end %ut
            if nout > 2, varargout(3)={-obj.a*cos(t) -   2*obj.b*cos(2*t)               };  end %utt
            if nout > 3, varargout(4)={ obj.a*sin(t) +   4*obj.b*sin(2*t)               };  end %u3t
            if nout > 4, varargout(5)={ obj.a*cos(t) +   8*obj.b*cos(2*t)               };  end %u4t

           
%             [Comp1{1:nout}] = obj.ACost.Derivatives( t );
%             [Comp2{1:nout}] = obj.HalfBCos2t.Derivatives(t);
%
%             varargout = cell(nout,1);
%            
%             for j=1:nout
%                 varargout(j) = { Comp1{j} + Comp2{j} };
%             end
%            
%             varargout(1) = {varargout{1} - obj.b/2};
        end
       
    end
end