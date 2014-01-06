classdef InteriorLaplacianHomoSolver < Solvers.SuperHomoSolver
   
    properties(Access = protected, AbortSet = true)
        Op;
        OpCoeffs;
    end
    
    methods
        function obj = InteriorLaplacianHomoSolver( ...
                BasisIndices,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams)
            obj = obj@Solvers.SuperHomoSolver( ...
                BasisIndices,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams);
            
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
%                ScattK = struct('r',abs(X+1i.*Y),'r0',CoeffsAddParams.r0);%ScattererClsHandle(GridK,obj.ScattererAddParams);
 %               WNPlr= Tools.Coeffs.CoeffsPolarR(ScattK,CoeffsAddParams.k0);
  %              obj.k = sparse(WNPlr.k);
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
        
    end
    
end