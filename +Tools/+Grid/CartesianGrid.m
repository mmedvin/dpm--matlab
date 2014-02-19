classdef CartesianGrid < Tools.Grid.Grids
    properties
        X;
        Y;        
        Z;
        R;
		Theta;
    end
    
    methods
        function obj = CartesianGrid(x1,xn,Nx,y1,yn,Ny)
            obj.x1 = x1;
            obj.xn = xn;
            obj.Nx = Nx;
            
            obj.y1 = y1;
            obj.yn = yn;
            obj.Ny = Ny;
            
            obj.dx=(xn-x1)/(Nx-1);
            obj.dy=(yn-y1)/(Ny-1);
            
            obj.x = obj.x1:obj.dx:obj.xn;
            obj.y = obj.y1:obj.dy:obj.yn;
        end
        
        function z = get.Z(obj)
            [X_,Y_] = obj.mesh();
            z = X_+1i*Y_;
        end
        
        function r = get.R(obj)
            r = abs(obj.Z());
		end
		
		function r = get.Theta(obj)
            r = angle(obj.Z());
        end
     
        function x = get.X(obj)
            x = meshgrid(obj.x,obj.y);
        end
        
        function y = get.Y(obj)
            [~,y] = meshgrid(obj.x,obj.y);
        end
    end
end
