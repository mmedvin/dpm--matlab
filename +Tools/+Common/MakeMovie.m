function MakeMovie(filename)
    
    Record = Input(filename);
    
    MovieName = filename;
    VideoWriterObj = VideoWriter([MovieName '.avi']);
    VideoWriterObj.FrameRate = 1/5;
    open(VideoWriterObj);
    
    FigNum=200;
    
    for indx=1:Record.FramesStored
        M=real(Record.Movie(indx).Mat);
        MAX(indx) = max(M(:));
        MIN(indx) = min(M(:));
    end
    
    MAX=max(MAX);
    MIN=min(MIN);
    
    for indx=1:1:Record.FramesStored
        
        figure(FigNum),
        axis([Record.Grid.x1,Record.Grid.xn,Record.Grid.y1,Record.Grid.yn, MIN,MAX])
        
        %%
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
        t=Record.Movie(indx).t;
        title(['t=' num2str(t)])
        %drawnow;
        
        writeVideo(VideoWriterObj,getframe(FigNum));
    end
    
    close(VideoWriterObj)
end


function Record = Input(filename)
    S = load(filename);

    Names = fieldnames(S);
    Record = eval(['S.' Names{1}]);
    
%     namesWorkspace = who;
%     outStr = regexpi(namesWorkspace, 'nameOfVariable');
%     ind = ~cellfun('isempty',outStr);
%     vars = namesWorkspace(ind);
end
