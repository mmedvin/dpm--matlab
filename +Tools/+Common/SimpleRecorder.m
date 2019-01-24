classdef SimpleRecorder < Tools.Common.AbstractRecorder
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FramesStored;
        Grid;
                
        Movie;
    end
    
    methods
        function obj = SimpleRecorder(Grid, NFrames)
            obj.FramesStored = 0;
            obj.Grid = Grid;
                
            Z = zeros(Grid.Nx,Grid.Ny);
            
            obj.Movie = struct('t',0,'Mat',Z);   % workaround a possible matlab bug, the following line won't work without this one as of matlab R2017a
            obj.Movie(1:NFrames) = struct('t',0,'Mat',Z);

        end
        
        function Store(obj, t, Mat)
            obj.FramesStored = obj.FramesStored + 1;
            
            obj.Movie(obj.FramesStored).t=t;
            obj.Movie(obj.FramesStored).Mat = Mat;
        end
    end
    
end

