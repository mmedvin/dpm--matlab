function varargout = FzeroFindMyZeros...
        (indx,InitialGuess,GridPointToTest,XHandle,YHandle)
    
    options = [];
    
    [varargout{1:nargout}] = fzero(@FMZ , GridPointToTest,options);
    
    
    
    function g = FindMyZeros(t,theta)
        
        x0 = XHandle.Derivatives(t);
        y0 = YHandle.Derivatives(t);
        
        arg = angle(x0+1i*y0);
        if sign(arg) < sign(t) && abs(theta) > pi/2
            arg = arg + 2*pi;
        elseif sign(arg) > sign(t) && abs(theta) > pi/2
            arg = arg - 2*pi;
        end
        
        g = arg - theta;
    end
    
    function varargout = FMZ(arg)
        [varargout{1:nargout}] = FindMyZeros(arg,InitialGuess);
    end
    
    
    
end






% function varargout = FzeroFindMyZeros...
%         (indx,InitialGuess,GridPointToTest,XHandle,YHandle)
%     
%     options = [];
%     
%     [varargout{1:nargout}] = fzero(@FMZ , GridPointToTest(indx),options);
%     
%     
%     
%     function g = FindMyZeros(t,theta)
%         
%         x0 = XHandle.Derivatives(t);
%         y0 = YHandle.Derivatives(t);
%         
%         arg = angle(x0+1i*y0);
%         if sign(arg) < sign(t) && abs(theta) > pi/2
%             arg = arg + 2*pi;
%         elseif sign(arg) > sign(t) && abs(theta) > pi/2
%             arg = arg - 2*pi;
%         end
%         
%         g = arg - theta;
%     end
%     
%     function varargout = FMZ(arg)
%         [varargout{1:nargout}] = FindMyZeros(arg,InitialGuess(indx));
%     end
%     
%     
%     
% end
