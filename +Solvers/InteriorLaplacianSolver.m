classdef InteriorLaplacianSolver < Solvers.SuperNonHomoSolver
   
    properties(Access = protected, AbortSet = true)
        Op;
        %OpCoeffs;
    end
    
    methods
        function obj = InteriorLaplacianSolver(Arguments)
            obj = obj@Solvers.SuperNonHomoSolver(Arguments);
            
%             if  numel(obj.Coeffs.a) + numel(obj.Coeffs.b) + numel(obj.Coeffs.sigma)>3
%                 GridK = Tools.Grid.CartesianGrid( ...
%                     obj.Grid.x1 - obj.Grid.dx/2 , ...
%                     obj.Grid.xn + obj.Grid.dx/2 , ...
%                     2*obj.Grid.Nx + 1           , ...
%                     obj.Grid.y1 - obj.Grid.dy/2 , ...
%                     obj.Grid.yn + obj.Grid.dy/2 , ...
%                     2*obj.Grid.Ny + 1         ) ;
%                                 
%                 ScattK = struct('r',GridK.R);
% 				obj.OpCoeffs=Tools.Coeffs.LaplaceCoeffsPolar(ScattK,CoeffsParams);
% 				
%             else
%                 obj.OpCoeffs = CoeffsParams;
%             end
                                
            %BCinRhs=0;
            %obj.Op =  Tools.DifferentialOps.LaplacianOp(Grid,obj.OpCoeffs,BCinRhs);
				
			Arguments.DiffOpParams.Grid             =   Arguments.Grid;
			Arguments.DiffOpParams.CoeffsHandle     =   Arguments.CoeffsHandle;
			Arguments.DiffOpParams.CoeffsParams     =   Arguments.CoeffsParams;
            if isfield(Arguments.ScattererParams,'FocalDistance')
                Arguments.DiffOpParams.CoeffsParams.FocalDistance = Arguments.ScattererParams.FocalDistance;
            end
            
            obj.Op = Arguments.DiffOp(Arguments.DiffOpParams);
        end
		
		function u = P_Omega(obj,xi_gamma)
			
			rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
			
			rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
			GLW = obj.Gf(rhs(:));
			
			u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));
			u(obj.Scatterer.Np)=xi_gamma(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).' + obj.GF{1}(obj.Scatterer.Np).';
		end
    end
    
    methods(Access = protected)
      
        function f = Lu(obj,u,msk)
            if isstruct(u), msk=u.msk; u = u.W; end
            
            if exist('msk','var');
				%f = obj.Op.A(msk,:)*u;
				f = obj.Op.ApplyOp(u,msk);
            else
				%f = obj.Op.A*u;
				f = obj.Op.ApplyOp(u);
             end
        end
        
        function u = Gf(obj,f)
			%u = obj.Op.A\f;
            if numel(f) == obj.Grid.Nx*obj.Grid.Ny
                u = obj.Op.Solve(f(:));
            else
			u = obj.Op.Solve(f);
        end
        end
                
        function u = GfSrc(obj,f)
            u = obj.Op.Solve(f(:),0);
        end

        
        function Qj = Qcol(obj,GLW,~)
            Qj = -GLW(obj.GridGamma,:);
		end
        
		function rhs = Bf(obj,F)
			%rhs = F(:);
			rhs = obj.Op.Bf(F);
		end
		
		function res = BF(obj)
			ScattererForSource = obj.Scatterer;
			
			HS = obj.SourceHandle(ScattererForSource,obj.CoeffsHandle,obj.CoeffsParams,obj.SourceParams);
			
			%res = obj.Bf(HS.Source);
            res = obj.Op.Bf(HS.Source);
			res(obj.Scatterer.Outside())=0;			
		end
	end
end
