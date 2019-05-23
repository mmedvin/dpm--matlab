classdef AbstractEnum
     methods
        function s= toString(obj)
            s= char(obj);
        end
        %         function s = char(obj)
        %             s = ['Color ' num2str(obj)];
        %             %# or use a switch statement..
        %         end
        %
        %         function disp(obj)
        %             disp( char(obj) )
        %         end
    end
end

