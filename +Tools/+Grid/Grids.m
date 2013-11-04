classdef Grids < handle
    
    properties(Access = public)
        x1;
        xn;
        Nx;
        dx;
        
        y1;
        yn;
        Ny;
        dy;
        
        x;% x-vector: x1:dx:xn
        y;% y-vector: y1:dy:yn
        
        Size;                
    end
    
    methods
        function obj = Grids(x1,xn,Nx,y1,yn,Ny)
            if nargin>0
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
        end
        
		function S = get.Size(obj) 
			S = [obj.Nx,obj.Ny];
		end

        function [X,Y] = mesh(obj)
            [X,Y] = meshgrid(obj.x,obj.y);
        end     
    end
    
end