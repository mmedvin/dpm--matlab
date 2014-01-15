classdef ExteriorLaplacianSolver < Solvers.InteriorLaplacianSolver

    methods
        function obj = ExteriorLaplacianSolver( ...
                Basis,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams,Source,SourceParams)
            obj = obj@Solvers.InteriorLaplacianSolver( ...
                Basis,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams,Source,SourceParams);            
        end
        
        function u = P_Omega(obj,xi_gamma)
            
            rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
            
            rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
            GLW = obj.Gf(rhs(:));
            
            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Nm));
			%u(obj.Scatterer.Nm)=xi_gamma(obj.Scatterer.Nm) - GLW(obj.Scatterer.Nm).';%( GLW(obj.Scatterer.Nm).' - xi_gamma(obj.Scatterer.Nm));
			u(obj.Scatterer.Nm)=xi_gamma(obj.Scatterer.Nm)  + ( GLW(obj.Scatterer.Nm).' - xi_gamma(obj.Scatterer.Nm)) + obj.GF(obj.Scatterer.Nm).';			
        end
    end
    
    methods(Access = protected)
                
        function Qj = Qcol(obj,GLW,w)
            Qj = GLW(obj.GridGamma,:)-w(obj.GridGamma,:);
        end
        

		
		function res = BF(obj)
			 
			ScattererForSource = obj.Scatterer;
			
			 HS = obj.SourceHandle(ScattererForSource,obj.CoeffsClsHandle,obj.CoeffsAddParams,obj.SourceParams);
			 
			 res = obj.Bf(HS.Source);
			 res(obj.Scatterer.Inside())=0;
			 
		 end
		
    end
        
end
