classdef SuperNonHomoSolver < Solvers.SuperHomoSolver
     properties(Access = public)
        Qf;
        Wf;
        GF;
        TrGF;
     end
     
     properties(Access = protected)
       %  Source;
         rhsf;  
         SourceHandle;
         myQf;
     end
     
     
     methods
         
          function qf = get.Qf(obj)
            if obj.IsReadyQnW == false
                obj.calc_QnW();
            end
            
            qf=obj.myQf;
            
        end
         
         function obj = SuperNonHomoSolver( ...
             Basis,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams,SourceHandle)
             
             obj = obj@Solvers.SuperHomoSolver( ...
                 Basis,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams);
             obj.SourceHandle = SourceHandle;
                                     
            % obj.Source = obj.SourceHandle(obj.Scatterer.TheScatterer,obj.WaveNumberClsHandle,obj.WaveNumberAddParams);
                          
             obj.rhsf=spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             obj.GF = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
             obj.Wf = spalloc(obj.Grid.Nx,obj.Grid.Ny,length(obj.GridGamma));

             
             
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
             
             obj.myQ0 = obj.Qcol( GLW(:,           1:obj.Basis.NBss   ) );
             obj.myQ1 = obj.Qcol( GLW(:,(obj.Basis.NBss+1):2*obj.Basis.NBss ) );
             
             
             obj.myQf = obj.Qcol( GLW(:,2*obj.Basis.NBss + 1 ) );
             
            obj.GF(obj.Scatterer.Np)  = GLW(obj.Scatterer.Np          ,2*obj.Basis.NBss + 2);
            obj.TrGF                = GLW(obj.Scatterer.GridGamma   ,2*obj.Basis.NBss + 2);
            
            obj.IsReadyQnW = true;
         end
         
         function Expand(obj)
             Expand@Solvers.SuperHomoSolver(obj);
             
             NoXi = obj.Basis.Handle();
             
           %  HS = obj.Source(obj.FocalDist,obj.Eta0,obj.phi,obj.k0,obj.r0); %?????
           
           Source = obj.SourceHandle(obj.Scatterer.TheScatterer,obj.WaveNumberClsHandle,obj.WaveNumberAddParams);
           
             obj.Wf(obj.GridGamma) = obj.Scatterer.Expansion(NoXi,NoXi,Source,obj.WaveNumber);             
         end
                  
        
                
        function res = BF(obj)
            %             global Eta x y x1 xn y1 yn dx dy cols rows IntEta k0 k
            
            % obj.WN = WaveNumberElliptical(FocalDist,eta,phi,k0,r0);
            % Exact = ExactElpsVarKr(FocalDist,eta, phi, obj.WN);
            
            
%             GridF = Grids( ...
%                         obj.Grid.x1 - obj.Grid.dx/2 , ...
%                         obj.Grid.xn + obj.Grid.dx/2 , ...
%                         obj.Grid.Nx * 2 + 1         , ...
%                         obj.Grid.y1 - obj.Grid.dy/2 , ...
%                         obj.Grid.yn + obj.Grid.dy/2 , ...
%                         obj.Grid.Ny * 2 + 1         ) ;
%             
%             
%             ScattF = obj.ScattererClsHandle(GridF,obj.ScattererAddParams);
           % WN = WaveNumberClsHandle(ScattF,obj.WaveNumberAddParams);
           % HS = SourceHandle(ScattF,WN);
           
           GridF = Tools.Grid.CartesianGrid( ...
                        obj.Grid.x1 - obj.Grid.dx/2 , ...
                        obj.Grid.xn + obj.Grid.dx/2 , ...
                        obj.Grid.Nx * 2 + 1         , ...
                        obj.Grid.y1 - obj.Grid.dy/2 , ...
                        obj.Grid.yn + obj.Grid.dy/2 , ...
                        obj.Grid.Ny * 2 + 1         ) ;
            ScattererForSource = obj.ScattererClsHandle(GridF,obj.ScattererAddParams);
           HS = obj.SourceHandle(ScattererForSource,obj.WaveNumberClsHandle,obj.WaveNumberAddParams);
            
            
            
% %             
% %             [X,Y] = meshgrid(obj.x(1)-obj.dx/2:obj.dx/2:obj.x(end)+obj.dx/2 ,obj.y(1)-obj.dy/2:obj.dy/2:obj.y(end)+obj.dy/2);
% %             Z=X+1i*Y;
% %             
% %             ETA = real(acosh(Z/obj.FocalDist));
% %             PHI = imag(acosh(Z/obj.FocalDist));
% %             
% %             GrdGamma = IntCartSlvrElps.SplitGrid(ETA, obj.Eta0);
%             
%             F=zeros(size(ETA));
%             %             F(GrdGamma)= approx_f(obj.FocalDist,obj.Eta0,PHI(GrdGamma),obj.k0,ETA(GrdGamma) - obj.Eta0);
%                         
%             HS = obj.Source(obj.FocalDist,obj.Eta0,PHI(GrdGamma),obj.k0,r0);
%             
%             %tmp   = calc_f(obj.FocalDist,ETA,PHI,obj.k0);
%             
%             DETA=ETA(GrdGamma) - obj.Eta0;
%             F(GrdGamma) = HS.F + DETA.*HS.Fn  + (DETA.^2).*HS.Fnn/2;%taylor
%             
%             tmp = obj.Source.calc_F(obj.FocalDist,ETA,PHI,obj.k0,r0);            
%             F(ETA<=obj.Eta0) = tmp(ETA<=obj.Eta0);
%             
%             %F=tmp;
%             BF  = IntCartSlvrElps.HlmSemBf(obj.x,obj.y,obj.k,F);
%             %obj.k=[];%kinda clean

           
            res = obj.Bf(HS.Source);
          %  BF=reshape(BF,obj.Grid.Nx,obj.Grid.Ny);
            res(obj.Scatterer.Outside())=0;
            
        end
        
        
        
        
     end
end
