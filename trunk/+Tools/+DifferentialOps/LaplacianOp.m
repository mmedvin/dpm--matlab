classdef LaplacianOp<Tools.DifferentialOps.SuperDiffOp
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
    properties(Access = public)
        A;
    end
    properties(Access = protected)
        Coeffs;
        Grid;
        BCinRhs;
		SolverType = 0;
		MG;
    end
    
    methods (Access=public)
        
        function obj = LaplacianOp(Grid, coeffs,BCinRhs,SolverType)
            obj.Grid=Grid;
            obj.VerifyCoeffs(coeffs);
            obj.BCinRhs=BCinRhs;
            obj.LaplacianA();
			
			switch SolverType
				case 1
					obj.SolverType = 1;
					
					mg_param.nu1 = 5;            % number of pre smooths
					mg_param.nu2 = 5;            % number of post smooths
					mg_param.minlevel = 1;       % level of coarsest grid
					mg_param.bx  = 5;            % block size x
					mg_param.by  = 5;            % block size y
					mg_param.sx  = 4;            % skip size x
					mg_param.sy  = 4;            % skip size y
					mg_param.cycleType = 'f';             % MG cycle type.  % the other choice is 'v'
					mg_param.verbose = 0;
					mg_param.maxIt = 100;
					
					if 1
						mg_param.TOL = (Grid.dx*Grid.dy)^2;
					else
						mg_param.TOL = 1.e-10;
					end
					
					nVoffset=-1;
					
					assert(Grid.Nx==Grid.Ny);
					obj.MG = Tools.LASolvers.MultiGrid.ClassQinghaiMG(obj.A, Grid.Nx, 'p1',mg_param, nVoffset);
				case 0
					% do nothing
				otherwise
					error('SolverType is not recognized')
			end
				
			end
				
		
		
			function res = Solve(obj,rhs)
				switch obj.SolverType
					case 0
						% spparms('spumoni',1)
						res = obj.A\rhs;
						%spparms('spumoni',0)
					case 1
						
						[res,resRel, nIters] = obj.MG.Solve(rhs);
				end
			end
			
        function Rhs = AdjustRhs(obj,Rhs,Exact)
                       
            Nx = obj.Grid.Nx;
            Ny = obj.Grid.Ny;
                
            if obj.BCinRhs                
				dx2 = obj.Grid.dx.^2;
				dy2 = obj.Grid.dy.^2;
                
                a_iphalf_j = obj.Coeffs.a(2:2:end-1,3:2:end);
                a_imhalf_j = obj.Coeffs.a(2:2:end-1,1:2:end-1);
                
                b_i_jphalf = obj.Coeffs.b(3:2:end,2:2:end-1);
                b_i_jmhalf = obj.Coeffs.b(1:2:end-1,2:2:end-1);                
                
                
                Rhs(1:Nx,   1   ) = Rhs(1:Nx, 1    ) - a_iphalf_j(:	, 1	).*Exact(2:end-1 , 1      )./dx2;%m
                Rhs(1:Nx,   Ny  ) = Rhs(1:Nx, Ny   ) - a_imhalf_j(:	, Ny).*Exact(2:end-1 , end    )./dx2;%p
                Rhs(1   ,  1:Ny ) = Rhs(1   , 1:Ny ) - b_i_jphalf(1	, : ).*Exact(1       , 2:end-1)./dy2; %m
                Rhs(Nx  ,  1:Ny ) = Rhs(Nx  , 1:Ny ) - b_i_jmhalf(Nx, :	).*Exact(end     , 2:end-1)./dy2; %p
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
                                
                a_iphalf_j = obj.Coeffs.a(2:2:end-1,3:2:end)./dx2;
                a_imhalf_j = obj.Coeffs.a(2:2:end-1,1:2:end-1)./dx2;
                
                b_i_jphalf = obj.Coeffs.b(3:2:end,2:2:end-1)./dy2;
                b_i_jmhalf = obj.Coeffs.b(1:2:end-1,2:2:end-1)./dy2;
                
                
                Aij = -( (a_iphalf_j(:) + a_imhalf_j(:))+ (b_i_jphalf(:) + b_i_jmhalf(:)) + obj.Coeffs.Sigma(:) );
                
                Aip1j = [zeros(1,cols)          ,a_iphalf_j(1:end-cols) ].';
                Aim1j = [a_imhalf_j(cols+1:end) ,zeros(1,cols)          ].';
                Aijp1 = [0                      ,b_i_jphalf(1:end-1)    ].';
                Aijm1 = [b_i_jmhalf(2:end)      ,0                      ].';
                
                Aijm1(cols:cols:rc)=0;
                Aijp1(cols+1:cols:rc)=0;
                
                obj.A = spdiags([Aim1j,Aijm1,Aij,Aijp1,Aip1j],[-cols,-1,0,1,cols], rc,rc);
                

            end
            
        end
          
        function VerifyCoeffs(obj,coeffs)
            if numel(coeffs.a)==1
                obj.Coeffs.a=  coeffs.a*ones( 2*obj.Grid.Nx + 1, 2*obj.Grid.Ny + 1 ) ;
            else
                obj.Coeffs.a = coeffs.a;
            end
            if numel(coeffs.b)==1
                obj.Coeffs.b=  coeffs.b*ones( 2*obj.Grid.Nx + 1, 2*obj.Grid.Ny + 1 ) ;
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
                                    
            if     any(size(obj.Coeffs.a) ~= [2*obj.Grid.Nx+1 , 2*obj.Grid.Ny+1]) || ...
                   any(size(obj.Coeffs.b) ~= [2*obj.Grid.Nx+1 , 2*obj.Grid.Ny+1])                   
                error('non constant coefficients a and b should be evaluated at x1-dx/2:dx:xn+dx/2, y1-dy/2:dy:yn+dy/2')
            end
            
            
%             obj.Coeffs.a = sparse(obj.Coeffs.a);
%             obj.Coeffs.b = sparse(obj.Coeffs.b);
%             obj.Coeffs.Sigma = sparse(obj.Coeffs.Sigma);
        end
        
    end
    
end

