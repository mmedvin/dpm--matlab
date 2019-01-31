classdef Lee2DRecorder < Tools.AbstractRecorder
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Omega;
        Psi;
        
        FramesStored;
                
        PrepareMovie
        Movie;
    end
    
    methods
        function obj  = NSRecorder(Grid, NFrames)
            obj.FramesStored = 0;
                
            Z = zeros(Grid.Nx,Grid.Ny);
            
            obj.Omega = struct('t',0,'Mat',Z);   % workaround a possible matlab bug, the following line won't work without this one as of matlab R2017a
            obj.Omega(1:NFrames) = struct('t',0,'Mat',Z);
            
            obj.Psi = struct('t',0,'Mat',Z);   % workaround a possible matlab bug, the following line won't work without this one as of matlab R2017a
            obj.Psi(1:NFrames) = struct('t',0,'Mat',Z); 
        end
        
        function         Store(obj, t, O,P)
            obj.FramesStored = obj.FramesStored + 1;
            
            obj.Omega(obj.FramesStored).t=t;
            obj.Omega(obj.FramesStored).Mat = O;
             
             obj.Psi(obj.FramesStored).t = t;
             obj.Psi(obj.FramesStored).Mat = P;
             
        end
    end
    
end

