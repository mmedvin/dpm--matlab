classdef PolarGrids < Tools.Grid.Grids
    properties
        r;
        theta;
        R;
        Theta;
        
        Z; %return X+1i*Z
    end
    methods
        function obj = PolarGrids(r1,rn,Nr,Nth)
            obj.x1 = r1;
            obj.xn = rn;
            obj.Nx = Nr;
            
            obj.dx=(rn-r1)/(Nr-1);
            obj.dy=2*pi/(Nth-1);%theta
            
            obj.y1 = 0;%+obj.dy/2;
            obj.yn = 2*pi - obj.dy;%/2;
            obj.Ny = Nth-1;
            
            obj.x = obj.x1:obj.dx:obj.xn;
            obj.y = obj.y1:obj.dy:obj.yn;
        end
        
        function z = get.Z(obj)
            [R_,Theta_] = obj.mesh();
            z = R_.*(cos(Theta_)+1i*sin(Theta_));
            z=z.';
        end
        
        function r_ = get.r(obj)
            r_ = obj.x;
        end
        
        function th_ = get.theta(obj)
            th_ = obj.y;
        end
        
        function r = get.R(obj)
            r = obj.mesh();
            r=r.';
        end
        
        function th = get.Theta(obj)
            [~,th] = obj.mesh();
            th=th.';
        end
        
    end
end
