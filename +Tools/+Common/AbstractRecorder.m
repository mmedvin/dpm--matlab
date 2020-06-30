classdef AbstractRecorder < handle
% abstract recorder class     
% recorder class supposed to save required pieces of data

    properties
    end
    
    methods(Abstract)
        Store(obj);%(obj, rho, u, v);
    end
    
end

