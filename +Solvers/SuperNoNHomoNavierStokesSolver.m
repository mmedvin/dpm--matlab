classdef SuperNoNHomoNavierStokesSolver < Solvers.SuperHomoNavierStokesSolver
     properties(SetAccess = protected)
		Qall;
		Wall;
        Qf;
        Wf;
        GF;
        TrGF;
        %GGF;
        TrGGF;
        QPsif;
        WPsif;
 
        TrGpsiGExf;
     end
     
     properties(Access = protected)
       %  Source;
         rhsf;  
         SourceHandle;
		 SourceParams=[];
         myQf;
         myQPsif

		 myGF;
		 myTrGF;
         myTrGPsiGF;
         rhsPsif;
         
         ExpansionPsiOrder;
     end
     
     methods(Access = protected,Abstract=true)
         %res = 
            Bf(obj,F);

         TrGPsiGf(obj,Gf);
     end
        
     methods 
         
         function UpdateSource(obj, Params)
             obj.SourceParams = Params.SourceParams;
                          
             obj.CreateWf();
             obj.calc_QnWf();
         end
         
		 function qf = get.Qf(obj)
			 if obj.IsReadyQnW == false
				 obj.calc_QnW();
			 end
			 
			 %qf=obj.myQf;
			 qf={obj.NewQ{:,end}};
         end
		 
         function qf = get.QPsif(obj)
             if obj.IsReadyQnW == false
                 obj.calc_QnW();
             end
             
             qf=obj.myQPsif;
         end

         function wf = get.WPsif(obj)
             if obj.IsReadyQnW == false
                 obj.calc_QnW();
             end
             
             wf=obj.ExtensionPsi.Wf;
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
             
             if isequal(Arguments.ExtensionPsi , @Tools.Extensions.NavierStokesPsi4rdOrderExtension)
                 obj.ExpansionPsiOrder =4;
             elseif isequal(obj.ExtensionPsi , @Tools.Extensions.NavierStokesPsi5rdOrderExtension)
                 obj.ExpansionPsiOrder=5;
             else
                 warning(' unknown ExtensionPsi!!!!');
             end
             
         end
     end
     
     methods(Access = protected)
         
         function Rhs(obj)
             Rhs@Solvers.SuperHomoNavierStokesSolver(obj);
             obj.CreateRhsf();
         end
         
         function calc_QnW(obj)
             
             calc_QnW@Solvers.SuperHomoNavierStokesSolver(obj);
             
             %obj.CreateWf();it is called inside SuperHomoNavierStokesSolver via abstract functions???
             [~,m2]=size(obj.Extension.Wf);
             
             obj.NewQ(1,end+(1:m2)) = {[]};
             obj.calc_QnWf();
             
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
             
             obj.calc_Qpsif();
         end
         
         function u = SolveSrc(obj,x)
             u = obj.Solve(x);
         end
         
         function Expand(obj)
             Expand@Solvers.SuperHomoNavierStokesSolver(obj);
             
             obj.CreateWf();                     
         end

         function calc_Qpsif(obj)

             if obj.ExpansionPsiOrder ==5
                 GLWPsif  = obj.Gf(obj.rhsPsif,obj.OpPsi);
                 
                 obj.myQPsif = obj.Qcol(GLWPsif,obj.ExtensionPsi.Wf);
             else
                 obj.myQPsif = spalloc(numel(obj.Scatterer.GridGamma),1,0);
             end             
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function CreateWf(obj)
           obj.Extension.ExpandSource(obj.SourceHandle,obj.SourceParams);
         
           if obj.ExpansionPsiOrder ==5
               obj.ExtensionPsi.ExpandSource(obj.SourceHandle,obj.SourceParams);
           end
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

            if obj.ExpansionPsiOrder ==5
                tmp= obj.Lu(obj.ExtensionPsi.Wf{1},obj.Scatterer.Mp,obj.OpPsi);
                obj.rhsPsif = obj.doitPsi(obj.ExtensionPsi.Wf{1},tmp);
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
