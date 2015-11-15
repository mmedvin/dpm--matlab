classdef Angles
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Radians;
        Degrees;
    end
    
    methods
        function  obj = Angles(value,type)
            if ~exist('type', 'var'), type=degrees; end
            
            switch type
                case 'radians'
                     obj. Radians = value;
                      obj. Degrees = value*180/pi;
                case 'degrees'
                      obj. Radians=value*pi/180;
                      obj. Degrees = value;
                otherwise 
                    error('unsupported degree type');
            end
        end
        
%         function rad = double(obj)
%             rad = obj.Radians;
%         end
%         
%         function disp(obj)
%             disp(obj.Degrees)
%         end
        
%         function str = char(obj)
%             str = sprintf('%s deg',obj.Degrees);
%         end
    end
    
end

