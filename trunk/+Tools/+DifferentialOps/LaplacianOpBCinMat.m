classdef LaplacianOpBCinMat<Tools.DifferentialOps.SuperLaplacianOp
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
    properties
        BC_x1=0;
        BC_xn=0;
        BC_y1=0;
        BC_yn=0;
    end
    methods (Access=public)
        
        function obj = LaplacianOpBCinMat(ParamsStruct)
            obj = obj@Tools.DifferentialOps.SuperLaplacianOp(ParamsStruct);	
			obj.ApplyBC();
            
            obj.BC_x1 = ParamsStruct.BC_x1;
            obj.BC_xn = ParamsStruct.BC_xn;
            obj.BC_y1 = ParamsStruct.BC_y1;
            obj.BC_yn = ParamsStruct.BC_yn;
            
        end
		function Rhs = AdjustRhs(obj,Rhs,Exact)
			if exist('Exact','var')
			Nx = obj.Grid.Nx;
			Ny = obj.Grid.Ny;
			
			Rhs(1:Nx,  1      ) = Exact(1:Nx,  1      );
			Rhs(1:Nx,  Ny     ) = Exact(1:Nx,  Ny     );
			Rhs(1   ,  1: Ny  ) = Exact(1   ,  1:Ny   );
			Rhs(Nx  ,  1: Ny  ) = Exact(Nx  ,  1:Ny   );
            else 
                
                Rhs(:	, 1		) = obj.BC_y1;%m
                Rhs(:   , end	) = obj.BC_yn;%p
                Rhs(1	, :     ) = obj.BC_x1; %m
                Rhs(end	, :     ) = obj.BC_xn; %p
			end
		end
    end
    methods (Access=protected)
        function ApplyBC(obj)
             
            rows = obj.Grid.Nx;
            cols = obj.Grid.Ny;
            
            rc=rows*cols;
   
			obj.A(1:cols,:)=  sparse(1:cols,1:cols,1,rows,rc);					% Ic; %1,j
			obj.A(1:cols:rc,:) = sparse(1:cols,1:cols:rc,1,rows,rc);			% i,1
			obj.A((rc-cols+1):rc,:)=sparse(1:cols,(rc-cols+1):rc,1,rows,rc);	% n,j
			obj.A(cols:cols:rc,:) = sparse(1:cols,cols:cols:rc,1,rows,rc);		%  i,n
   
        end
    end
    
end

