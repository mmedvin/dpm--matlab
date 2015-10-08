classdef SuperNonHomoSolver < Solvers.SuperHomoSolver
     properties(Access = public)
		Qall;
		Wall;
        Qf;
        Wf;
        GF;
        TrGF;
     end
     
     properties(Access = protected)
       %  Source;
         rhsf;  
         SourceHandle;
		 SourceParams=[];
         myQf;

		 myGF;
		 myTrGF;
     end
     
     
     methods
         
		 function qf = get.Qf(obj)
			 if obj.IsReadyQnW == false
				 obj.calc_QnW();
			 end
			 
			 %qf=obj.myQf;
			 qf=obj.NewQ{end};
		 end
		 
		 function wf = get.Wf(obj)
			 if obj.IsReadyQnW == false
				 obj.calc_QnW();
			 end
			 
			 wf=obj.Extension.Wf;
		 end
		 
		 function res = get.TrGF(obj)
			 if obj.IsReadyQnW == false
				 obj.calc_QnW();
			 end
			 
			 res=obj.myTrGF;
			 
		 end
		 
		 function res = get.GF(obj)
			 if obj.IsReadyQnW == false
				 obj.calc_QnW();
			 end
			 
			 res=obj.myGF;
			 
		 end
         		 
         function obj = SuperNonHomoSolver(Arguments)
             
             obj = obj@Solvers.SuperHomoSolver(Arguments);

             obj.SourceHandle = Arguments.SourceHandle;
             if isfield(Arguments,'SourceParams')
                 obj.SourceParams = Arguments.SourceParams;
             end
             
             obj.rhsf=spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             obj.myGF = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             
         end
     end
     
     methods(Abstract, Access = protected)
        %res = 
            Bf(obj,F);
     end
     
     methods(Access = protected)
         
         function Rhs(obj)
             Rhs@Solvers.SuperHomoSolver(obj);
             obj.CreateRhsf();
         end
         
         function calc_QnW(obj)
			 if obj.CollectRhs
				 obj.Rhs();
				 				 
                 GLW = cellfun(@(arg) obj.Gf(arg),[obj.rhs,{obj.rhsf(:)},{obj.BF}],'UniformOutput',false);
                 obj.NewQ = cellfun(@(arg1,arg2) obj.Qcol(arg1,arg2),GLW(1:end-1), [obj.Extension.W,{obj.Extension.Wf}],'UniformOutput',false);
                 
                 obj.myGF   = GLW{end};
                 obj.myTrGF = GLW{end}(obj.Scatterer.GridGamma);
			 else
                 calc_QnW@Solvers.SuperHomoSolver(obj);
                 
                 obj.CreateWf();
                 
                 GLW = obj.SolveSrc(obj.Extensions.Wf);                 
                 obj.myQf = obj.Qcol(GLW,obj.Extensions.Wf);
                 
                 obj.myGF   = obj.Gf(obj.BF);
                 obj.myTrGF = obj.myGF(obj.Scatterer.GridGamma);
			 end
			 obj.IsReadyQnW = true;
         end
         
         function u = SolveSrc(obj,x)
             u = obj.Solve(x);
%              obj.f(obj.Scatterer.Mp) = obj.Lu(x(:),obj.Scatterer.Mp);
%              %obj.Truncate(f);
%              u = obj.Gf(obj.f);%(:));
         end
         
         function Expand(obj)
             Expand@Solvers.SuperHomoSolver(obj);
             
             obj.CreateWf();                     
         end
                  
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function CreateWf(obj)
             %NoXi = obj.Basis.Handle();
             
           %  HS = obj.Source(obj.FocalDist,obj.Eta0,obj.phi,obj.k0,obj.r0); %?????
           
           Source = obj.SourceHandle(obj.Scatterer.TheScatterer,obj.CoeffsHandle,obj.CoeffsParams,obj.SourceParams);
           
           %obj.myWf(obj.GridGamma) = obj.Scatterer.Expansion(NoXi,NoXi,Source,obj.Coeffs);
           obj.Extension.ExpandSource(Source);
         end
                  
         function CreateRhsf(obj)             
             obj.rhsf(obj.Scatterer.Mp) = obj.Lu(obj.Extension.Wf(:),obj.Scatterer.Mp);
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
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
			 HS = obj.SourceHandle(ScattererForSource,obj.CoeffsHandle,obj.CoeffsParams,SourceParamsForSource);
			 
			 res = obj.Bf(HS.Source);
			 res(obj.Scatterer.Outside())=0;
			 
		 end
        
        
        
        
     end
end
