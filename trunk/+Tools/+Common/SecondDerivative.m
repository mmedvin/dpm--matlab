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
            
            Ym2 = sparse(3:obj.cols,1:obj.cols-2,1,obj.cols,obj.cols);
            Yp2=Ym2';
            
            %Y = Ym1+Yp1;
            
            Xm1 = sparse(2:obj.rows,1:obj.rows-1,1,obj.rows,obj.rows);
            Xp1=Xm1';
            %X = Xm1+Xp1;
            
            Xm2 = sparse(3:obj.rows,1:obj.rows-2,1,obj.rows,obj.rows);
            Xp2=Xm2';
            
            Ay = kron( Ir                   , (Yp1 - Ym1)/2/obj.dy    );
            Ax = kron( (Xp1 - Xm1)/2/obj.dx , Ic                    );

            ur  = (Ax*u ).*obj.dxdr    + (Ay*u).*obj.dydr;
            
            Ayy = kron( Ir                     , (Yp1 + Ym1)/(obj.dy^2)    ) -2*II/(obj.dy^2);
            Axx = kron( (Xp1 + Xm1)/(obj.dx^2) , Ic                    ) -2*II/(obj.dx^2);
            %%Axy = kron( (Xp1 - Xm1)/(obj.dx*obj.dy) , Ic ) +  kron( Ir , (Yp1 - Ym1)/(obj.dx*obj.dy)    );
            %%Axy = - kron( (Xp1 + Xm1)/obj.dx , (Yp1 + Ym1)/obj.dy ) +  kron( Ir , (Yp1 - Ym1)/(obj.dx*obj.dy)    );
          
            
            %Ay = kron( Ir                                       , (-Yp2 + 8*Yp1 - 8*Ym1 + Ym2)/12/obj.dy    );
            %Ax = kron( (-Xp2 + 8*Xp1 - 8*Xm1 + Xm2)/12/obj.dx   , Ic                                        );
            
            %Ayy = kron( Ir                          , (-Yp2 + 16*Yp1 + 16*Ym1 - Ym2)/12/(obj.dy^2)     ) -30*II/12/(obj.dy^2);
            %Axx = kron( (-Xp2 + 16*Xp1 + 16*Xm1 - Xm2)/12/(obj.dx^2)   , Ic                            ) -30*II/12/(obj.dx^2);
            
            urr = (Axx*u).*obj.dxdr.^2 + (Ayy*u).*obj.dydr.^2 + 2*((Ax*Ay)*u).*obj.dxdr.*obj.dydr;
            
        end
        
    end
    
end