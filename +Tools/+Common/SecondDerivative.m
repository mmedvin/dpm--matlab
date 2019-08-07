classdef SecondDerivative < handle
    
    properties
        rows;
        cols;
        dx;
        dy;
        dxdr=1;
        dydr=1;
        
        Ax;
        Ay;
        Axx;
        Ayy;
        Axy;
    end
    
    
    methods
        
        function obj = SecondDerivative(rows,cols,dx,dy,dxdr,dydr)
            obj.rows = rows;
            obj.cols = cols;
            obj.dx   = dx;
            obj.dy   = dy;
            
            if nargin > 4
                obj.dxdr = dxdr;
                obj.dydr = dydr;
            end
            
            obj.CreateOps();
        end
            
        function CreateOps(obj)
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
            
            obj.Ay = kron( Ir                   , (Yp1 - Ym1)/2/obj.dy    );
            obj.Ax = kron( (Xp1 - Xm1)/2/obj.dx , Ic                    );

            %%%%%%%%%%%%%%%%%
            
            obj.Ayy = kron( Ir                     , (Yp1 + Ym1)/(obj.dy^2)    ) -2*II/(obj.dy^2);
            obj.Axx = kron( (Xp1 + Xm1)/(obj.dx^2) , Ic                    ) -2*II/(obj.dx^2);
            obj.Axy = obj.Ax*obj.Ay;
          
            
            %%%%%%%%%%%%%%%% Fourth Order %%%%%%%%%%%%%%%%
            %Ay = kron( Ir                                       , (-Yp2 + 8*Yp1 - 8*Ym1 + Ym2)/12/obj.dy    );
            %Ax = kron( (-Xp2 + 8*Xp1 - 8*Xm1 + Xm2)/12/obj.dx   , Ic                                        );
            
            %Ayy = kron( Ir                          , (-Yp2 + 16*Yp1 + 16*Ym1 - Ym2)/12/(obj.dy^2)     ) -30*II/12/(obj.dy^2);
            %Axx = kron( (-Xp2 + 16*Xp1 + 16*Xm1 - Xm2)/12/(obj.dx^2)   , Ic                            ) -30*II/12/(obj.dx^2);
            
            

        end
        
        function [ux,uy,uxx,uyy,uxy] = CartesianDerivatives(obj,u)

            ux = obj.Ax*u(:);
            uy = obj.Ay*u(:);
            uxx= obj.Axx*u(:);
            uyy= obj.Ayy*u(:);
            uxy= obj.Axy*u(:);

            ux  = reshape(ux  , obj.rows , obj.cols);
            uy  = reshape(uy  , obj.rows , obj.cols);
            uxx = reshape(uxx , obj.rows , obj.cols);
            uyy = reshape(uyy , obj.rows , obj.cols);
            uxy = reshape(uxy , obj.rows , obj.cols);

        end
        
        function L = CartesianLaplacian(obj,u)
            [~,~,uxx,uyy] = obj.CartesianDerivatives(u);
            L = uxx + uyy;
        end
         
        function [ur,urr] = RadialDerivatives(obj,u)            

            ur  = (obj.Ax*u(:) ).*obj.dxdr    + (obj.Ay*u(:)).*obj.dydr;
                        
            urr = (obj.Axx*u(:)).*obj.dxdr.^2 + (obj.Ayy*u(:)).*obj.dydr.^2 + 2*(obj.Axy*u(:)).*obj.dxdr.*obj.dydr;
            
            ur  = reshape(ur  , obj.rows , obj.cols);
            urr = reshape(urr , obj.rows , obj.cols);
        end
        
    end
    
end