function MakeMovie(filename)
    
    load(filename)
    
    MovieName = filename;
    VideoWriterObj = VideoWriter([MovieName '.avi']);
    VideoWriterObj.FrameRate = 1/5;
    open(VideoWriterObj);
    
    FigNum=200;
    
    for indx=1:Record.FramesStored
        
        figure(FigNum),
        M=real(Record.Movie(indx).Mat);
        
        mesh(Record.Grid.X,Record.Grid.Y,M);
        %             if nargin==2
        %                  pcolor(Grid.X,Grid.Y,M);
        %             else
        %                 %imagesc(Grid.X,Grid.Y,M);
        %                 imagesc(M);
        %             end
        
        
        %-->
        colorbar
        %axis equal
        title(['t=' num2str(Record.Movie(indx).t)])
        %drawnow;
        
        writeVideo(VideoWriterObj,getframe(FigNum));
    end
    
    close(VideoWriterObj)
end

