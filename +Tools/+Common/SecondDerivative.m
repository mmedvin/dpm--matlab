classdef SecondDerivative < handle
    
    properties
        rows;
        cols;
        dx;
        dy;
        dxdr;
        dydr;
    end
    
    
    methods
        
        function obj = SecondDerivative(rows,cols,dx,dy,dxdr,dydr)
            obj.rows = rows;
            obj.cols = cols;
            obj.dxdr = dxdr;
            obj.dydr = dydr;
            obj.dx   = dx;
            obj.dy   = dy;
        end
        function [ur,urr] = RadialDerivatives(obj,u)
            
            Ir = speye(obj.rows);
            Ic = speye(obj.cols);
            II = speye(obj.rows*obj.cols);
            
            Ym1 = sparse(2:obj.cols,1:obj.cols-1,1,obj.cols,obj.cols);
            Yp1=Ym1';
            %Y = Ym1+Yp1;
            
            Xm1 = sparse(2:obj.rows,1:obj.rows-1,1,obj.rows,obj.rows);
            Xp1=Xm1';
            %X = Xm1+Xp1;
            
            Ay = kron( Ir                   , (Yp1 - Ym1)/2/obj.dy    );
            Ax = kron( (Xp1 - Xm1)/2/obj.dx , Ic                    );

            Ayy = kron( Ir                     , (Yp1 + Ym1)/(obj.dy^2)    ) -2*II/(obj.dy^2);
            Axx = kron( (Xp1 + Xm1)/(obj.dx^2) , Ic                    ) -2*II/(obj.dx^2);
            %Axy = kron( (Xp1 - Xm1)/(obj.dx*obj.dy) , Ic ) +  kron( Ir , (Yp1 - Ym1)/(obj.dx*obj.dy)    );
            %Axy = - kron( (Xp1 + Xm1)/obj.dx , (Yp1 + Ym1)/obj.dy ) +  kron( Ir , (Yp1 - Ym1)/(obj.dx*obj.dy)    );
          
            
            ur  = (Ax*u ).*obj.dxdr    + (Ay*u).*obj.dydr;
            urr = (Axx*u).*obj.dxdr.^2 + (Ayy*u).*obj.dydr.^2 + 2*((Ax*Ay)*u).*obj.dxdr.*obj.dydr;
            
        end
        
    end
    
end