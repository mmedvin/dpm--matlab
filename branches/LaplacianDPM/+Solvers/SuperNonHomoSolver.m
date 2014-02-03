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
		 SourceParams;
         myQf;
     end
     
     
     methods
         
		 function qf = get.Qf(obj)
			 if obj.IsReadyQnW == false
				 obj.calc_QnW();
			 end
			 
			 qf=obj.myQf;
			 
		 end
         
		 function q = get.Qall(obj)
            q=[obj.Q,obj.myQf];
		 end
		 
		 function w = get.Wall(obj)
			 w=[obj.W,obj.Wf(:)];
		 end
		 
         function obj = SuperNonHomoSolver( ...
             Basis,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams,SourceHandle,SourceParams)
             
             obj = obj@Solvers.SuperHomoSolver( ...
                 Basis,Grid,CoeffsClsHandle,CoeffsAddParams,ScattererClsHandle,ScattererAddParams);
             obj.SourceHandle = SourceHandle;
			  if exist('SourceParams','var')
				  obj.SourceParams = SourceParams;
			  else
				  obj.SourceParams=[];
			  end
			 
                                     
            % obj.Source = obj.SourceHandle(obj.Scatterer.TheScatterer,obj.WaveNumberClsHandle,obj.WaveNumberAddParams);
                          
             obj.rhsf=spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             obj.GF = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             obj.Wf = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.GridGamma));

             
             
         end
     end
     
     methods(Abstract, Access = protected)
        %res = 
            Bf(obj,F);
     end
     
     methods(Access = protected)
         
         function Rhs(obj)
             Rhs@Solvers.SuperHomoSolver(obj);
             tmp =obj.Lu(obj.Wf(:));
             obj.rhsf(obj.Scatterer.Mp) = tmp(obj.Scatterer.Mp,:);
         end
         
         function calc_QnW(obj)
             obj.Rhs();
                          
             GLW = obj.Gf([obj.rhs0,obj.rhs1,obj.rhsf(:),obj.BF]);
             
             obj.myQ0 = obj.Qcol( GLW(:,           1:obj.Basis.NBss   )		 , obj.W0 );
             obj.myQ1 = obj.Qcol( GLW(:,(obj.Basis.NBss+1):2*obj.Basis.NBss ), obj.W1 );
			 
			 obj.myQf = obj.Qcol( GLW(:,2*obj.Basis.NBss + 1 )				 , obj.Wf(:) );
             
            %obj.GF(obj.Scatterer.Np)  = GLW(obj.Scatterer.Np          ,2*obj.Basis.NBss + 2);
			obj.GF  = GLW(:         ,2*obj.Basis.NBss + 2);
            obj.TrGF                = GLW(obj.Scatterer.GridGamma   ,2*obj.Basis.NBss + 2);
            
            obj.IsReadyQnW = true;
         end
         
         function Expand(obj)
             Expand@Solvers.SuperHomoSolver(obj);
             
             NoXi = obj.Basis.Handle();
             
           %  HS = obj.Source(obj.FocalDist,obj.Eta0,obj.phi,obj.k0,obj.r0); %?????
           
           Source = obj.SourceHandle(obj.Scatterer.TheScatterer,obj.CoeffsClsHandle,obj.CoeffsAddParams,obj.SourceParams);
           
             obj.Wf(obj.GridGamma) = obj.Scatterer.Expansion(NoXi,NoXi,Source,obj.Coeffs);             
         end
                  
        
                
		 function res = BF(obj)
			 
			 GridF = Tools.Grid.CartesianGrid( ...
				 obj.Grid.x1 - obj.Grid.dx/2 , ...
				 obj.Grid.xn + obj.Grid.dx/2 , ...
				 obj.Grid.Nx * 2 + 1         , ...
				 obj.Grid.y1 - obj.Grid.dy/2 , ...
				 obj.Grid.yn + obj.Grid.dy/2 , ...
				 obj.Grid.Ny * 2 + 1         ) ;
			 ScattererForSource = obj.ScattererClsHandle(GridF,obj.ScattererAddParams);
			 HS = obj.SourceHandle(ScattererForSource,obj.CoeffsClsHandle,obj.CoeffsAddParams);
			 
			 res = obj.Bf(HS.Source);
			 res(obj.Scatterer.Outside())=0;
			 
		 end
        
        
        
        
     end
end
