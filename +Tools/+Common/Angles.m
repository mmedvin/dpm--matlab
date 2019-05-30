classdef Angles
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        Radians;
        Degrees;
    end
    
        properties(Access = protected)
            ValueInDegrees;
        end
    
    methods
        function  obj = Angles(value,type)
           % if ~exist('type', 'var'), type=Tools.Enums.Angles.Degrees; end

            switch type
                case Tools.Enums.Angles.Radians
                    obj.ValueInDegrees = value*180/pi;
                case Tools.Enums.Angles.Degrees
                    obj. ValueInDegrees = value;
                otherwise
                    error('unsupported degree type');
            end
        end
        
        
        
        function rad = get.Radians(obj)
            rad = obj.ValueInDegrees*pi/180;
        end
        
        function deg = get.Degrees(obj)
            deg = obj.ValueInDegrees;
        end
        
%         function val=double(obj)
%             val = obj.Radians;
%         end
%         
%         function disp(obj)
%             disp(obj.Degrees)
%         end
%         
%         function str = char(obj)
%             str = sprintf('%d deg',obj.ValueInDegrees);
%         end
    end
    
end

