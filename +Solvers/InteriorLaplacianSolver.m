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
                    obj.Grid.x1 - obj.Grid.dx/2 , ...
                    obj.Grid.xn + obj.Grid.dx/2 , ...
                    2*obj.Grid.Nx + 1           , ...
                    obj.Grid.y1 - obj.Grid.dy/2 , ...
                    obj.Grid.yn + obj.Grid.dy/2 , ...
                    2*obj.Grid.Ny + 1         ) ;
                                
                ScattK = struct('r',GridK.R);
				obj.OpCoeffs=Tools.Coeffs.LaplaceCoeffsPolar(ScattK,CoeffsAddParams);
				
            else
                obj.OpCoeffs = CoeffsAddParams;
            end
            
            BCinRhs=0;
            obj.Op =  Tools.DifferentialOps.LaplacianOp(Grid,obj.OpCoeffs,BCinRhs);
        end
		
		function u = P_Omega(obj,xi_gamma)
			
			rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
			
			rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
			GLW = obj.Gf(rhs(:));
			
			u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));
			u(obj.Scatterer.Np)=xi_gamma(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).' + obj.GF(obj.Scatterer.Np).';
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
			rhs = F(:);
			% do nothing, required for compact scheme only 
		end
		
		function res = BF(obj)
			
			ScattererForSource = obj.Scatterer;
			
			HS = obj.SourceHandle(ScattererForSource,obj.CoeffsClsHandle,obj.CoeffsAddParams,obj.SourceParams);
			
			res = obj.Bf(HS.Source);
% 			Src  = HS.Source;
% 			res = Src;
 			res(obj.Scatterer.Outside())=0;
% 			
% 			res(1:end,1)	= Src(1:end,1);
%             res(1,1:end)	= Src(1,1:end);
%             res(1:end,end)	= Src(1:end,end);
%             res(end,1:end)	= Src(end,1:end);
% 			res = res(:);
			
		end
		
    end
        
end
