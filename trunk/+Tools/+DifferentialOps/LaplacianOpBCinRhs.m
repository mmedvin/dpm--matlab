classdef LaplacianOpBCinRhs<Tools.DifferentialOps.SuperLaplacianOp
    %LaplacianOp Creates 2nd order Matrix for variable coefficiant Laplacian
   
	properties (Constant)
		GigaByte=2^30;
	end
	
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
		L;
		U;
		P;
		Q;
		R;
		tollerance = 1e-10; %default value, may be changed in some iterative solvers
		MaxIteration ;
        RestartAt;
        
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
						mg_param.minlevel = 4;       % level of coarsest grid
						mg_param.bx  = 5;            % block size x
						mg_param.by  = 5;            % block size y
						mg_param.sx  = 4;            % skip size x
						mg_param.sy  = 4;            % skip size y
						mg_param.cycleType = 'v';    % MG cycle type.  % the other choice is 'v', 'f'
						mg_param.verbose = 0;
						mg_param.maxIt = 100;
						
						if 1
							mg_param.TOL = max((EGrid.dx*EGrid.dy)^2,obj.tollerance);
						else
							mg_param.TOL = obj.tollerance; %1.e-10;
						end
						
						nVoffset=-1;
						
						assert(ParamsStruct.Grid.Nx==ParamsStruct.Grid.Ny);
						obj.MGHandle = Tools.LASolvers.MultiGrid.ClassQinghaiMG(obj.A, ParamsStruct.Grid.Nx, 'p1',mg_param, nVoffset);
				case 2 %MATAMG, aka AMG_MP
                    
                    obj.tollerance =  (obj.Grid.dx*obj.Grid.dy);%^2;
					
					Params.MaxCycle = 100;
					Params.tollerance = obj.tollerance;
                    Params.A = obj.A;
					 obj.MGHandle = Tools.LASolvers.MultiGrid.MatAmgWrapper(Params);
					
                    %obj.AMG_MP_Options = amgset;
                    %obj.AMG_MP_Options = amgset(obj.AMG_MP_Options,'PrintOnScreen','off', 'MaxCycle',200,'PreCond','pcg','TolAMG', obj.tollerance);%'IntpType','lramg');%, 'SaveCsn','off','SaveIntp','off','Log','off');
					
				case 3 %LU
					[obj.L,obj.U,obj.P,obj.Q,obj.R] = lu(obj.A);					
				case 4 %LU	
					[obj.L,obj.U,obj.P,obj.Q] = lu(obj.A);
					
				case 5 % GMRES
					
					obj.tollerance		=  (obj.Grid.dx*obj.Grid.dy)^2; %max((obj.Grid.dx*obj.Grid.dy)^2,1e-8);
					obj.MaxIteration	=	obj.Grid.Nx*obj.Grid.Ny;
					
					
					%[obj.L,obj.U] = ilu(obj.A,struct('type','nofill'));%preconditioner
					[obj.L,obj.U] = ilu(obj.A,struct('type','ilutp','droptol',obj.tollerance));%preconditioner
					%[obj.L,obj.U] = ilu(obj.A,struct('type','crout','droptol',obj.tollerance));%preconditioner
					
                    tmp = obj.A;
                    s=whos('tmp');
                    
                    obj.RestartAt = 2*fix(obj.GigaByte/s.bytes); % the system should have at least 2gbyte, we want to restart on a half way
                    maxit = floor(obj.MaxIteration/obj.RestartAt);
                    if  obj.RestartAt<1, obj.RestartAt=10; maxit = floor(obj.MaxIteration/obj.RestartAt); end
                    if obj.RestartAt> obj.Grid.Nx*obj.Grid.Ny, obj.RestartAt = []; maxit = obj.MaxIteration; end %means no restart
                    
                    obj.MaxIteration = maxit;
                
                case 6
                    
                    obj.tollerance		=  (obj.Grid.dx*obj.Grid.dy)^2; %max((obj.Grid.dx*obj.Grid.dy)^2,1e-8);
					obj.MaxIteration	=	sqrt(obj.Grid.Nx*obj.Grid.Ny);
					
					
					%[obj.L,obj.U] = ilu(obj.A,struct('type','nofill'));%preconditioner
					[obj.L,obj.U] = ilu(obj.A,struct('type','ilutp','droptol',obj.tollerance));%preconditioner
					%[obj.L,obj.U] = ilu(obj.A,struct('type','crout','droptol',obj.tollerance));%preconditioner

				case 99 %AMG_WU
					
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
		
		function u = Solve(obj,f,SolverType)
			u = zeros(size(f));
			Rhs = f(obj.Inside);
			
            if ~exist('SolverType','var')
                SolverType = obj.LinearSolverType;
            end
            
			%for j=1:size(f,2)
				switch SolverType 
					case 1 % Multigrid of Qinghai
						%assert(size(f,2)==1);
						%[u(obj.Inside,j),resRel, nIters] = obj.MGHandle.Solve(f(obj.Inside,j));
						[u(obj.Inside),resRel, nIters] = obj.MGHandle.Solve(Rhs);
											
                    case 2 %MATAMG, aka AMG_MP
                        
						u(obj.Inside) = obj.MGHandle.Solve(Rhs);                
						
                      % u(obj.Inside)= amg(obj.A,u(obj.Inside),f(obj.Inside),obj.AMG_MP_Options);
                                              
					  %                        b=f(obj.Inside);
					  %                        B=obj.A*u(obj.Inside);
					  %                        assert(norm(b(:)-B(:),'inf')<1e-10);

					case 3
						u(obj.Inside)=obj.Q*(obj.U\(obj.L\(obj.P*(obj.R\Rhs))));
					case 4
						u(obj.Inside)=obj.Q*(obj.U\(obj.L\(obj.P*(Rhs))));

					case 5 % GMRES
						%assert(size(f,2)==1);
						tic
                        
						InitialGuess = obj.U\(obj.L\Rhs);
						%InitialGuess = rand(size(f(obj.Inside)));
						
                        [u(obj.Inside), flag,relres,iter,resvec] = gmres(obj.A ,f(obj.Inside), obj.RestartAt, obj.tollerance,obj.MaxIteration,obj.L,obj.U,InitialGuess);
												
						t1=toc; 
						if flag, %warning('gmres return with flag=%d \n',flag); 	
							fprintf('gmres done with flag=%d,outer iter=%d, inner iter=%d, tol=%d, relres =%d time=%d \n',flag,iter(1),iter(2), tol, relres, t1);
						end
                        
                    case 6
                        InitialGuess = obj.U\(obj.L\Rhs);
                        
                        [u(obj.Inside),flag,relres] = bicgstab(obj.A,f(obj.Inside),obj.tollerance,obj.MaxIteration,obj.L,obj.U,InitialGuess);
                        
					case 99 %AMG_WU
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

