classdef LaplacianOpBCinRhs<Tools.DifferentialOps.SuperLaplacianOp
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
    properties
		Inside;
		nx;
		ny;
		
		BC_x1=0;
		BC_xn=0;
		BC_y1=0;
		BC_yn=0;
	end
    methods (Access=public)
        
        function obj = LaplacianOpBCinRhs(ParamsStruct)
			
			EGrid = ParamsStruct.Grid;			
			ParamsStruct.Grid = Tools.Grid.CartesianGrid(	EGrid.x(1)+EGrid.dx,EGrid.x(end)-EGrid.dx,EGrid.Nx-2, ...
															EGrid.y(1)+EGrid.dy,EGrid.y(end)-EGrid.dy,EGrid.Ny-2);
			
            obj = obj@Tools.DifferentialOps.SuperLaplacianOp(ParamsStruct);	
			
			obj.nx = EGrid.Nx;
			obj.ny = EGrid.Ny;
			rc = obj.nx*obj.ny;
			
			[I,J]=meshgrid(2:obj.nx-1,2:obj.ny-1);
			obj.Inside = sub2ind(size(EGrid.X),J,I);
			obj.Inside = obj.Inside(:);
			
			obj.BC_x1 = ParamsStruct.BC_x1;
			obj.BC_xn = ParamsStruct.BC_xn;
			obj.BC_y1 = ParamsStruct.BC_y1;
			obj.BC_yn = ParamsStruct.BC_yn;
		end
		
				function b = ApplyOp(obj,x,mask)
			% returns matrix-vector multiplication result, i.e b = A*x
			% if parameter 'mask' is esixts, a mask 'mask' is applied on the returned value b,
			% more precisely b(mask) is returned, i.e. only the indices in 'mask'
			
			b = zeros(size(x));
			b(obj.Inside,:) = obj.A*x(obj.Inside,:);
			if exist('mask','var');
				b=b(mask,:);
			end
		end
		
		function u = Solve(obj,f)
			u = zeros(size(f));
			
			for j=1:size(f,2)
				u(obj.Inside,j) = obj.A\f(obj.Inside,j);
			end
        end
		
		
		
        function Rhs = AdjustRhs(obj,Rhs,Exact,BC_yn,BC_x1,BC_xn)
                       
            Nx = obj.Grid.Nx;
            Ny = obj.Grid.Ny;
			
			dx2 = obj.Grid.dx.^2;
			dy2 = obj.Grid.dy.^2;
			
			a_iphalf_j = obj.Coeffs.a(2:2:end-1,3);
			a_imhalf_j = obj.Coeffs.a(2:2:end-1,end-1);
			
			b_i_jphalf = obj.Coeffs.b(1,2:2:end-1);
			b_i_jmhalf = obj.Coeffs.b(end-1,2:2:end-1);
			
			
% 			if nargin == 3
% 				BC_y1 = Exact(2:end-1 , 1      );
% 				BC_yn = Exact(2:end-1 , end    );
% 				BC_x1 = Exact(1       , 2:end-1);
% 				BC_xn = Exact(end     , 2:end-1);
% 				
% 			elseif nargin == 6
% 				BC_y1 = Exact;
% 			end
			
			Rhs(2:end-1	, 2			) = Rhs(2:end-1	, 2			) - a_iphalf_j.* obj.BC_y1./dx2;%m
			Rhs(2:end-1	, end-1		) = Rhs(2:end-1	, end-1		) - a_imhalf_j.* obj.BC_yn./dx2;%p
			Rhs(2		,  2:end-1	) = Rhs(2		, 2:end-1	) - b_i_jphalf.* obj.BC_x1./dy2; %m
			Rhs(end-1	,  2:end-1	) = Rhs(end-1	, 2:end-1	) - b_i_jmhalf.* obj.BC_xn./dy2; %p

			
			%old
			% 			a_iphalf_j = obj.Coeffs.a(2:2:end-1,3:2:end);
			% 			a_imhalf_j = obj.Coeffs.a(2:2:end-1,1:2:end-1);
			%
			% 			b_i_jphalf = obj.Coeffs.b(3:2:end,2:2:end-1);
			% 			b_i_jmhalf = obj.Coeffs.b(1:2:end-1,2:2:end-1);
			%
			%
			% 			Rhs(1:Nx,   1   ) = Rhs(1:Nx, 1    ) - a_iphalf_j(:	, 1	).*Exact(2:end-1 , 1      )./dx2;%m
			% 			Rhs(1:Nx,   Ny  ) = Rhs(1:Nx, Ny   ) - a_imhalf_j(:	, Ny).*Exact(2:end-1 , end    )./dx2;%p
			% 			Rhs(1   ,  1:Ny ) = Rhs(1   , 1:Ny ) - b_i_jphalf(1	, : ).*Exact(1       , 2:end-1)./dy2; %m
			% 			Rhs(Nx  ,  1:Ny ) = Rhs(Nx  , 1:Ny ) - b_i_jmhalf(Nx, :	).*Exact(end     , 2:end-1)./dy2; %p

        end
    end
 
    
end

