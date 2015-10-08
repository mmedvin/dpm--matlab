classdef ExteriorHomoLaplacianSolver < Solvers.InteriorHomoLaplacianSolver

	properties
% 		BC_y1;
% 		BC_yn;
% 		BC_x1;
% 		BC_xn;
	end
	
    methods
        function obj = ExteriorHomoLaplacianSolver(Arguments)
            obj = obj@Solvers.InteriorHomoLaplacianSolver(Arguments);     			
        end
        
        function u = P_Omega(obj,xi_gamma)
            
            rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
            
            rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
            GLW = obj.Gf(rhs(:));
            
            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Nm));
            %u(obj.Scatterer.Nm)=xi_gamma(obj.Scatterer.Nm)  + ( GLW(obj.Scatterer.Nm).' - xi_gamma(obj.Scatterer.Nm)) ;
            u(obj.Scatterer.Nm)=GLW(obj.Scatterer.Nm).' ;
        end
    end
    
    methods(Access = protected)
                
        function Qj = Qcol(obj,GLW,w)
            Qj = GLW(obj.GridGamma,:)-w(obj.GridGamma,:);
        end
        
	end
end
