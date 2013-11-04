classdef Angles
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Radians;
        Degrees;
    end
    
    methods
        function  obj = Angles(value,type)
            switch type
                case 'radians'
                     obj. Radians = value;
                      obj. Degrees = value*180/pi;
                case 'degrees'
                      obj. Radians=value*pi/180;
                      obj. Degrees = value;
            end
        end
    end
    
end

