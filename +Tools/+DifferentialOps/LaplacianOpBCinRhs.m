classdef LaplacianOpBCinRhs<Tools.DifferentialOps.SuperLaplacianOp
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
    properties
		
		LinearSolverType = 0;
		
		Inside;
		nx;
		ny;
		
		BC_x1=0;
		BC_xn=0;
		BC_y1=0;
		BC_yn=0;
		
		
		MGHandle;
		iluL;
		iluU;
		tollerance = 1e-12; %default value, may be changed in some iterative solvers
		gmresRestarts;
        
        AMG_MP_Options;
	end
    methods (Access=public)
        
        function obj = LaplacianOpBCinRhs(ParamsStruct)
			
			EGrid = ParamsStruct.Grid;			
			ParamsStruct.Grid = Tools.Grid.CartesianGrid(	EGrid.x(1)+EGrid.dx,EGrid.x(end)-EGrid.dx,EGrid.Nx-2, ...
															EGrid.y(1)+EGrid.dy,EGrid.y(end)-EGrid.dy,EGrid.Ny-2);
			
            obj = obj@Tools.DifferentialOps.SuperLaplacianOp(ParamsStruct);	
			
			obj.nx = EGrid.Nx;
			obj.ny = EGrid.Ny;
			
			
			[I,J]=meshgrid(2:obj.nx-1,2:obj.ny-1);
			obj.Inside = sub2ind(size(EGrid.X),J,I);
			obj.Inside = obj.Inside(:);
			
			obj.BC_x1 = ParamsStruct.BC_x1;
			obj.BC_xn = ParamsStruct.BC_xn;
			obj.BC_y1 = ParamsStruct.BC_y1;
			obj.BC_yn = ParamsStruct.BC_yn;
			
			
			obj.LinearSolverType = ParamsStruct.LinearSolverType;
			
			switch obj.LinearSolverType
					case 1 % Multigrid of Qinghai
						mg_param.nu1 = 3;            % number of pre smooths
						mg_param.nu2 = 3;            % number of post smooths
						mg_param.minlevel = 1;       % level of coarsest grid
						mg_param.bx  = 5;            % block size x
						mg_param.by  = 5;            % block size y
						mg_param.sx  = 4;            % skip size x
						mg_param.sy  = 4;            % skip size y
						mg_param.cycleType = 'v';    % MG cycle type.  % the other choice is 'v', 'f'
						mg_param.verbose = 0;
						mg_param.maxIt = 100;
						
						if 1
							mg_param.TOL = (EGrid.dx*EGrid.dy)^2;
						else
							mg_param.TOL = obj.tollerance; %1.e-10;
						end
						
						nVoffset=-1;
						
						assert(ParamsStruct.Grid.Nx==ParamsStruct.Grid.Ny);
						obj.MGHandle = Tools.LASolvers.MultiGrid.ClassQinghaiMG(obj.A, ParamsStruct.Grid.Nx, 'p1',mg_param, nVoffset);
					case 2 % GMRES
						
						obj.tollerance =  (obj.Grid.dx*obj.Grid.dy)^2;
						obj.gmresRestarts = 2*fix(sqrt(obj.Grid.Nx*obj.Grid.Ny));
						
						%[obj.iluL,obj.iluU] = ilu(obj.A,struct('type','nofill','droptol',obj.tollerance));%preconditioner
						[obj.iluL,obj.iluU] = ilu(obj.A,struct('type','ilutp','droptol',obj.tollerance^2));%preconditioner
						%[obj.iluL,obj.iluU] = ilu(obj.A,struct('type','crout','droptol',obj.tollerance));%preconditioner
										
                case 3 %MATAMG, aka AMG_MP
                    
                    obj.tollerance =  (obj.Grid.dx*obj.Grid.dy);%^2;
					
					Params.MaxCycle = 100;
					Params.tollerance = obj.tollerance;
                    Params.A = obj.A;
					 obj.MGHandle = Tools.LASolvers.MultiGrid.MatAmgWrapper(Params);
					
                    %obj.AMG_MP_Options = amgset;
                    %obj.AMG_MP_Options = amgset(obj.AMG_MP_Options,'PrintOnScreen','off', 'MaxCycle',200,'PreCond','pcg','TolAMG', obj.tollerance);%'IntpType','lramg');%, 'SaveCsn','off','SaveIntp','off','Log','off');
					
				
				case 4 %AMG_WU
					
					error('wrong choice, this option doesn''t work well')
					
					obj.tollerance =  (obj.Grid.dx*obj.Grid.dy)^2;
					
					Params.level=2;
					Params.relax_it=1;%2;
					Params.relaxation_parameter=1;%1
					Params.post_smoothing=1;
					Params.max_iter=100;
					Params.tollerance=obj.tollerance;
					Params.pc_type=2;
					Params.connection_threshold=0.25;
					
					obj.MGHandle = Tools.LASolvers.MultiGrid.AMG_WU(obj.A,Params);
				
				otherwise
						%nothing
				end
			
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
			
			%for j=1:size(f,2)
				switch obj.LinearSolverType
					case 1 % Multigrid of Qinghai
						assert(size(f,2)==1);
						%[u(obj.Inside,j),resRel, nIters] = obj.MGHandle.Solve(f(obj.Inside,j));
						[u(obj.Inside),resRel, nIters] = obj.MGHandle.Solve(f(obj.Inside));
					case 2 % GMRES
						assert(size(f,2)==1);
						[u(obj.Inside), flag] = gmres(obj.A ,f(obj.Inside), [], obj.tollerance, 100,obj.iluL,obj.iluU);
						
                    case 3 %MATAMG, aka AMG_MP
                        
						u(obj.Inside) = obj.MGHandle.Solve(f(obj.Inside));                
						
                      % u(obj.Inside)= amg(obj.A,u(obj.Inside),f(obj.Inside),obj.AMG_MP_Options);
                                              
%                        b=f(obj.Inside);
%                        B=obj.A*u(obj.Inside);
%                        assert(norm(b(:)-B(:),'inf')<1e-10);

					case 4 %AMG_WU
						assert(size(f,2)==1);
						[u(obj.Inside),resd, nIters] = obj.MGHandle.Solve(f(obj.Inside));                        

					otherwise
						%u(obj.Inside,j) = obj.A\f(obj.Inside,j);
						u(obj.Inside,:) = obj.A\f(obj.Inside,:);
				end
			%end
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

