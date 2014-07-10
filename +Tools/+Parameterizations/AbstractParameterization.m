classdef AbstractParameterization < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Access = public)
        XHandle;
        YHandle;
    end        
    
    methods
        function str = Print(obj)
            
            X = obj.XHandle;
            Y = obj.YHandle;
            Xmeta = metaclass(X);
            Ymeta = metaclass(Y);
                 
            str = sprintf('XHandle:%s, params: ',Xmeta.Name);
            
            for i = 1:numel(Xmeta.Properties)
                name = Xmeta.Properties{i}.Name;
                val = eval(['obj.XHandle.',name]);
                if isobject(val), continue,end
                
                str = sprintf('%s %s=%d',str,name,val);
            end
            
            str = sprintf('%s \n YHandle:%s, params: ',str, Ymeta.Name);
            
            for j = 1:numel(Ymeta.Properties)
                                
                name = Ymeta.Properties{j}.Name;
                val =  eval(['obj.YHandle.',name]);
                if isobject(val), continue,end
               
                str = sprintf('%s %s=%d',str,name,val);
            end

            %Tools.Parameterizations.
        end
        
%         function CompareTo(AnotherParametrization)
%             
%             obj.XHandle.CompareTo(AnotherParametrization.XHandle);
%             obj.YHandle.CompareTo(AnotherParametrization.YHandle);
%         end
    end
end

