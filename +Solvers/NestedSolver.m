classdef NestedSolver < Solvers.InteriorSolver
   
    properties(Access = protected, AbortSet = true)

    end
    
    methods
        function obj = NestedSolver(Arguments)
            obj = obj@Solvers.InteriorSolver(Arguments);                        
        end        
    end
    
    methods(Access = protected)
		 function res = BF(obj)
			 
			 GridF = Tools.Grid.CartesianGrid( ...
				 obj.Grid.x1 - obj.Grid.dx/2 , ...
				 obj.Grid.xn + obj.Grid.dx/2 , ...
				 obj.Grid.Nx * 2 + 1         , ...
				 obj.Grid.y1 - obj.Grid.dy/2 , ...
				 obj.Grid.yn + obj.Grid.dy/2 , ...
				 obj.Grid.Ny * 2 + 1         ) ;
			 ScattererForSource = obj.ScattererHandle(GridF,obj.ScattererParams);
             SourceParamsForSource = obj.SourceParams;
             SourceParamsForSource.Grid = GridF;
			 HSInt = obj.SourceHandle(ScattererForSource.InteriorScatterer,obj.CoeffsHandle,obj.CoeffsParams,SourceParamsForSource);
			 HSExt = obj.SourceHandle(ScattererForSource.ExteriorScatterer,obj.CoeffsHandle,obj.CoeffsParams,SourceParamsForSource);
			 
%              SrcInt= obj.Bf(HSInt.Source);
%              SrcInt(obj.Scatterer.InteriorScatterer.Outside()) = 0;
%              
%              SrcExt = obj.Bf(HSExt.Source);
%              SrcExt(obj.Scatterer.ExteriorScatterer.Outside()) = 0;

            Src = HSExt.Source-HSInt.Source;
            Src(ScattererForSource.InteriorScatterer.GridGamma) = HSInt.Source(ScattererForSource.InteriorScatterer.GridGamma);
            Src= obj.Bf(Src);
            Src(obj.Scatterer.Outside())=0;
            
			 res = {Src;Src };
			 %res(obj.Scatterer.Outside())=0;
			 
		 end        
    end
    
    
end
