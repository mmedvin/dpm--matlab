classdef InteriorLaplacianSolver < Solvers.SuperNonHomoSolver
   
    properties(Access = protected, AbortSet = true)
        Op;
        OpCoeffs;
    end
    
    methods
        function obj = InteriorLaplacianSolver( ...
                Basis,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams,Source,SourceParams)
            obj = obj@Solvers.SuperNonHomoSolver( ...
                Basis,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams,Source,SourceParams);
            
            if  numel(obj.Coeffs.a) + numel(obj.Coeffs.b) + numel(obj.Coeffs.sigma)>3
                GridK = Tools.Grid.CartesianGrid( ...
                    obj.Grid.x1 - obj.Grid.dx , ...
                    obj.Grid.xn + obj.Grid.dx , ...
                    obj.Grid.Nx + 2           , ...
                    obj.Grid.y1 - obj.Grid.dy , ...
                    obj.Grid.yn + obj.Grid.dy , ...
                    obj.Grid.Ny + 2         ) ;
                
                [X,Y] = GridK.mesh();
                
                error('TBD')
                %  ScattK = struct('r',abs(X+1i.*Y),'r0',CoeffsAddParams.r0);%ScattererClsHandle(GridK,obj.ScattererAddParams);
                %  WNPlr= Tools.Coeffs.CoeffsPolarR(ScattK,CoeffsAddParams.k0);
                %  obj.k = sparse(WNPlr.k);
            else
                obj.OpCoeffs = CoeffsAddParams;
            end
            
            BCinRhs=0;
            obj.Op =  Tools.DifferentialOps.LaplacianOp(Grid,obj.Coeffs,BCinRhs);
        end
        
        function u = P_Omega(obj,xi_gamma)
            
            rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
            
            rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
            GLW = obj.Gf(rhs(:));
            
            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
            u(obj.Scatterer.Np)=xi_gamma(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).';                        
        end
    end
    
    methods(Access = protected)
      
        function f = Lu(obj,u,msk)
            if exist('msk','var');
                f = obj.Op.A(msk,:)*u;
            else
                f = obj.Op.A*u;
            end
        end
        
        function u = Gf(obj,f)
           u = obj.Op.A\f;           
        end
                
        function Qj = Qcol(obj,GLW,~)
            Qj = -GLW(obj.GridGamma,:);
		end
        
		function rhs = Bf(obj,F)
			%obj.Op.AdjustRhs
% 			 F = F.';% probably no need
			rhs = F(:);
			% do nothing, required for compact scheme only 
		end
		
		function res = BF(obj)
			 
			% 			 GridF = Tools.Grid.CartesianGrid( ...
			% 				 obj.Grid.x1 - obj.Grid.dx/2 , ...
			% 				 obj.Grid.xn + obj.Grid.dx/2 , ...
			% 				 obj.Grid.Nx * 2 + 1         , ...
			% 				 obj.Grid.y1 - obj.Grid.dy/2 , ...
			% 				 obj.Grid.yn + obj.Grid.dy/2 , ...
			% 				 obj.Grid.Ny * 2 + 1         ) ;
			% 			 ScattererForSource = obj.ScattererClsHandle(GridF,obj.ScattererAddParams);
			ScattererForSource = obj.Scatterer;
			
			 HS = obj.SourceHandle(ScattererForSource,obj.CoeffsClsHandle,obj.CoeffsAddParams,obj.SourceParams);
			 
			 res = obj.Bf(HS.Source);
% 			 res(obj.Scatterer.Outside())=0;
			% res(obj.Scatterer.Inside())=0;
			 
		 end
		
    end
        
end
