classdef SuperNoNHomoNavierStokesSolver < Solvers.SuperHomoNavierStokesSolver
     properties(Access = public)
		Qall;
		Wall;
        Qf;
        Wf;
        GF;
        TrGF;
        %GGF;
        TrGGF;
        Qpsi_f;
        Wpsi_f; %???
        TrGpsiGExf;
     end
     
     properties(Access = protected)
       %  Source;
         rhsf;  
         SourceHandle;
		 SourceParams=[];
         myQf;

		 myGF;
		 myTrGF;
         myTrGPsiGF;
     end
     
     methods(Access = protected,Abstract=true)
         %res = 
            Bf(obj,F);

         TrGPsiGf(obj,Gf);
     end
        
     methods      
		 function qf = get.Qf(obj)
			 if obj.IsReadyQnW == false
				 obj.calc_QnW();
			 end
			 
			 %qf=obj.myQf;
			 qf={obj.NewQ{:,end}};
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
         	
         function res = get.TrGGF(obj)
             if obj.IsReadyQnW == false
                 obj.calc_QnW();
             end
             
             res=obj.myTrGPsiGF;
             
         end
         
         function obj = SuperNoNHomoNavierStokesSolver(Arguments)
             
             obj = obj@Solvers.SuperHomoNavierStokesSolver(Arguments);

             obj.SourceHandle = Arguments.SourceHandle;
             if isfield(Arguments,'SourceParams')
                 obj.SourceParams = Arguments.SourceParams;
             end
             
             obj.rhsf = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             obj.myGF = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             
         end
     end
     
     methods(Access = protected)
         
         function Rhs(obj)
             Rhs@Solvers.SuperHomoNavierStokesSolver(obj);
             obj.CreateRhsf();
         end
         
         function calc_QnW(obj)
% 			 if obj.CollectRhs
% 				 obj.Rhs();
% 				 				 
%                  GLW = cellfun(@(arg) obj.Gf(arg),[obj.rhs,obj.rhsf,obj.BF],'UniformOutput',false);
%                  obj.NewQ = cellfun(@(arg1,arg2) obj.Qcol(arg1,arg2),GLW(:,1:end-1), [obj.Extension.W,obj.Extension.Wf],'UniformOutput',false);
%                  
%                  obj.mP   = cellfun(@(arg1,arg2) obj.Pcol(arg1,arg2),GLW(:,1:end-2), obj.Extension.W,'UniformOutput',false);
%                  
%                  obj.myGF   = GLW(:,end);
%                   obj.myTrGPsiGF = obj.TrGPsiGf(obj.myGF{1});
%                   
%                  for indx=1:numel(GLW(:,end))
%                      obj.myTrGF{indx} = GLW{indx,end}(obj.Scatterer.GridGamma);
%                  end
% 			 else
                 calc_QnW@Solvers.SuperHomoNavierStokesSolver(obj);
                 
                 obj.CreateWf();
                 [~,m2]=size(obj.Extension.Wf);
%                  for i=1:m
%                     obj.NewQ{end+1} = obj.Extension.Wf(i);
%                  end
                obj.NewQ(1,end+(1:m2)) = {[]};%repmat({[]},m1);%perhaps numel should be put here instead of size????
                 obj.calc_QnWf();
                 
% 			 end
% 			 obj.IsReadyQnW = true;
         end         
         
         function calc_QnWf(obj)
             
             obj.myGF   = {obj.Gf(obj.BF)};
             obj.myTrGPsiGF = obj.TrGPsiGf(obj.myGF{1});
             
             m=numel(obj.Extension.Wf);
             for indx=1:m
                 GLW = obj.SolveSrc(obj.Extension.Wf{indx});
                 obj.NewQ{(end-m)+indx} = obj.Qcol(GLW,obj.Extension.Wf{indx});
                 obj.myTrGF{indx} = obj.myGF{1}(obj.Scatterer.GridGamma);
                 b = obj.TrGpsiPOmega(GLW,obj.Extension.Wf{indx});
                 obj.TrGpsiGExf = b;%(obj.GridGamma);
             end
             
             %m=numel(obj.NewQ);
             %  for indx=1:numel(obj.Extension.Wf)
             %      GLW = obj.SolveSrc(obj.Extension.Wf{indx});
             %      obj.NewQ{m+indx} = obj.Qcol(GLW,obj.Extension.Wf{indx});
             %      obj.myTrGF{indx} = obj.myGF{1}(obj.Scatterer.GridGamma);
             %  end
 
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
           
          % Source = obj.SourceHandle(obj.Scatterer.TheScatterer,obj.CoeffsHandle,obj.CoeffsParams,obj.SourceParams);
           
           %obj.myWf(obj.GridGamma) = obj.Scatterer.Expansion(NoXi,NoXi,Source,obj.Coeffs);
           %obj.Extension.ExpandSource(Source);
           obj.Extension.ExpandSource(obj.SourceHandle,obj.SourceParams);
           %obj.ExtensionPsi.ExpandSource(obj.SourceHandle,obj.SourceParams);
         end
                  
         function CreateRhsf(obj)             
             %obj.rhsf(obj.Scatterer.Mp) = obj.Lu(obj.Extension.Wf(:),obj.Scatterer.Mp);
             
             tmp=cellfun(@(arg) obj.Lu(arg,obj.Scatterer.Mp),obj.Extension.Wf,'UniformOutput',false);
            
            obj.rhsf = cell(size(tmp));
            for indx=1:numel(tmp)
                [n,m]=size(obj.Extension.Wf{indx});
                NNZ = nnz(obj.Extension.Wf{indx});
                obj.rhsf{indx} = spalloc( n,m,NNZ);
                obj.rhsf{indx}(obj.Scatterer.Mp,:) = tmp{indx};
            end

             
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
		 function Src = BF(obj)
			 
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
             Src={res};
			 
		 end
        
        
        
        
     end
end
