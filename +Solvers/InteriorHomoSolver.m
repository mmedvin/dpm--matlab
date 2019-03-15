classdef InteriorHomoSolver < Solvers.SuperHomoSolver
   
    properties(Access = protected, AbortSet = true)
       % A;
       % k;
       Op;
    end
    
    methods
        function obj = InteriorHomoSolver(Arguments)
            obj = obj@Solvers.SuperHomoSolver(Arguments);
            
            Arguments.DiffOpParams.Grid             =   Arguments.Grid;
            Arguments.DiffOpParams.CoeffsHandle     =   Arguments.CoeffsHandle;
            Arguments.DiffOpParams.CoeffsParams     =   Arguments.CoeffsParams;
   
            obj.Op = Arguments.DiffOp(Arguments.DiffOpParams);

            
        end
        
        function u = P_Omega(obj,xi_gamma)
            
            rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
            
            rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
            GLW = obj.Gf(rhs(:));
            
            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));
            u(obj.Scatterer.Np)=xi_gamma(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).';                        
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

                
        function Qj = Qcol(obj,GLW,W)
            Qj = -GLW(obj.GridGamma,:);
        end
        
           
    end
    
    
end
