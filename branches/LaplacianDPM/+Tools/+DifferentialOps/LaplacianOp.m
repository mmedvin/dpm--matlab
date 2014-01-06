classdef LaplacianOp<Tools.DifferentialOps.SuperDiffOp
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
    properties(Access = public)
        A;
    end
    properties(Access = protected)
        Coeffs;
        Grid;
        BCinRhs;
    end
    
    methods (Access=public)
        
        function obj = LaplacianOp(Grid, coeffs,BCinRhs)
            obj.Grid=Grid;
            obj.VerifyCoeffs(coeffs);
            obj.BCinRhs=BCinRhs;
            obj.LaplacianA();
        end
        function Rhs = AdjustRhs(obj,Rhs,Exact)
                       
            Nx = obj.Grid.Nx;
            Ny = obj.Grid.Ny;
                
            if obj.BCinRhs                
				dx2 = obj.Grid.dx.^2;
				dy2 = obj.Grid.dy.^2;
                Rhs(1:Nx,   1   ) = Rhs(1:Nx, 1    ) - obj.Coeffs.a(1:end-2 , 1       ).*Exact(2:end-1 , 1      )./dx2;%m
                Rhs(1:Nx,   Ny  ) = Rhs(1:Nx, Ny   ) - obj.Coeffs.a(3:end   , Ny+1    ).*Exact(2:end-1 , end    )./dx2;%p
                Rhs(1   ,  1:Ny ) = Rhs(1   , 1:Ny ) - obj.Coeffs.b(1       , 1:end-2 ).*Exact(1       , 2:end-1)./dy2; %m
                Rhs(Nx  ,  1:Ny ) = Rhs(Nx  , 1:Ny ) - obj.Coeffs.b(Nx      , 3:end   ).*Exact(end     , 2:end-1)./dy2; %p
            else
                Rhs(1:Nx,  1      ) = Exact(1:Nx,  1      );
                Rhs(1:Nx,  Ny     ) = Exact(1:Nx,  Ny     );
                Rhs(1   ,  1: Ny  ) = Exact(1   ,  1:Ny   );
                Rhs(Nx  ,  1: Ny  ) = Exact(Nx  ,  1:Ny   );                
            end
        end
    end
    methods (Access=protected)
        function LaplacianA(obj)%(x,y,k)
            dx = obj.Grid.dx;
            dy = obj.Grid.dy;
             
            rows = obj.Grid.Nx;
            cols = obj.Grid.Ny;
            
            rc=rows*cols;
                        
            InitialDerivationMatrixA();
                       
            %external boundary
            if obj.BCinRhs
%                 obj.A(1:cols,cols+1:end)=  spalloc(cols,rc-cols,0); %1,j
%                 obj.A(1:cols:rc,rc-cols-1:end) = spalloc(cols,rc-cols,0); %i,1
%                 obj.A((rc-cols+1):rc,1:rc-cols)= spalloc(cols,rc-cols,0); % n,j
%                 obj.A(cols:cols:rc,rc-cols:end) = spalloc(cols,rc-cols-1,0);%sparse(1:cols,cols:cols:rc,1,rows,rc); %  i,n
            else
                obj.A(1:cols,:)=  sparse(1:cols,1:cols,1,rows,rc);%Ic; %1,j
                obj.A(1:cols:rc,:) = sparse(1:cols,1:cols:rc,1,rows,rc); %i,1
                obj.A((rc-cols+1):rc,:)=sparse(1:cols,(rc-cols+1):rc,1,rows,rc); % n,j
                obj.A(cols:cols:rc,:) = sparse(1:cols,cols:cols:rc,1,rows,rc); %  i,n
            end
            
            function InitialDerivationMatrixA()
                % This function return  matrix A - inner scheme for Helmholtz equation
                
                dx2 = dx.^2;
                dy2 = dy.^2;
                                
                a_iphalf_j=obj.Coeffs.a(2:end-1,3:end);
                b_i_jphalf=obj.Coeffs.b(2:end-1,3:end);
                a_imhalf_j=obj.Coeffs.a(2:end-1,1:end-2);
                b_i_jmhalf=obj.Coeffs.b(2:end-1,1:end-2);
                
                
                Aij = -( (a_iphalf_j(:) + a_imhalf_j(:))./dx2+ (b_i_jphalf(:) + b_i_jmhalf(:))./dy2 + obj.Coeffs.Sigma(:) );
                
                Aim1j = a_imhalf_j(:)./dx2;
                Aim1j(cols:cols:rc) = 0;
                
                Aip1j = a_iphalf_j(:)./dx2;
                Aip1j(cols+1:cols:rc) = 0;
                
                Aijm1 = b_i_jphalf(:)./dy2;
                Aijp1 = b_i_jmhalf(:)./dy2;
                
                obj.A = spdiags([Aijm1, Aim1j,Aij,Aip1j,Aijp1],[-cols,-1,0,1,cols], rc,rc);
            end
            
        end
          
        function VerifyCoeffs(obj,coeffs)
            if numel(coeffs.a)==1
                obj.Coeffs.a=  coeffs.a*ones( obj.Grid.Nx + 2, obj.Grid.Ny + 2 ) ;
            else
                obj.Coeffs.a = coeffs.a;
            end
            if numel(coeffs.b)==1
                obj.Coeffs.b=  coeffs.b*ones( obj.Grid.Nx + 2, obj.Grid.Ny + 2 ) ;
            else
                obj.Coeffs.b = coeffs.b;
            end
            
            if numel(coeffs.sigma)==1
                obj.Coeffs.Sigma=  coeffs.sigma*ones( obj.Grid.Nx, obj.Grid.Ny) ;
            else 
                obj.Coeffs.Sigma  = coeffs.sigma;
            end
            
            
            if any(size(obj.Coeffs.Sigma) ~= [obj.Grid.Nx, obj.Grid.Ny])
                error('non constant sigma should be evaluated at x1:xn, y1:yn')
            end
                                    
            if     any(size(obj.Coeffs.a) ~= [obj.Grid.Nx+2 , obj.Grid.Ny+2]) || ...
                   any(size(obj.Coeffs.b) ~= [obj.Grid.Nx+2,obj.Grid.Ny+2])                   
                error('non constant coefficients a and b should be evaluated at x1-dx/2:xn+dx/2, y1-dy/2:yn+dy/2')
            end
        end
        
    end
    
end

