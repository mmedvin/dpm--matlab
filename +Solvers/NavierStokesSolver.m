classdef NavierStokesSolver < Solvers.SuperNoNHomoNavierStokesSolver
    %NAVIERSTOCKSHOMOSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = protected, AbortSet = true)
        Op;
        
        xi0Psi = 0;
        xi1Psi = 0;
        xi0PsiTT = 0;
    end
    properties(Access = public)%???
        OpPsi;
    end
    
      methods
          
          function UpdateSource(obj, Params)
             % obj.SourceHandle = NewSource;
              obj.SourceParams = Params.SourceParams;
              %               obj.Extension.ExpandSource(obj.SourceHandle,obj.SourceParams);
              %
              %               obj.CreateRhsf();
              %
              %               obj.CreateWf();
              %               obj.myGF   = {obj.Gf(obj.BF)};
              %               obj.myTrGPsiGF = obj.TrGPsiGf(obj.myGF{1});
              
              obj.xi0Psi  = Params.PsiBC.xi0Psi;
              obj.xi1Psi  = Params.PsiBC.xi1Psi;
              obj.xi0PsiTT= Params.PsiBC.xi0PsiTT;
              
              obj.CreateWf();
              obj.calc_QnWf();
              
              
 

          end
                
        function obj = NavierStokesSolver(Arguments)
            %SUPERHOMONAVIERSTOCKSSOLVER Construct an instance of this class
            %   Detailed explanation goes here
            
            obj = obj@Solvers.SuperNoNHomoNavierStokesSolver(Arguments);
            
            Arguments.DiffOpParams.Grid         = Arguments.Grid;
            Arguments.DiffOpParams.CoeffsHandle = Arguments.CoeffsHandle;
            Arguments.DiffOpParams.CoeffsParams = Arguments.CoeffsParams;
            
            if isfield(Arguments.ScattererParams,'FocalDistance')
                Arguments.DiffOpParams.CoeffsParams.FocalDistance = Arguments.ScattererParams.FocalDistance;
            end
            
            obj.Op = Arguments.DiffOp(Arguments.DiffOpParams);
            
            ArgPsi = Arguments.DiffOpParams;
            ArgPsi.CoeffsParams.sigma=0;
            %obj.OpPsi = Tools.DifferentialOps.LaplacianOpBCinRhs(ArgPsi);
            %Tools.DifferentialOps.LaplacianOp4(ArgPsi);
            %obj.OpPsi = Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs(ArgPsi);
            obj.OpPsi = Arguments.DiffOp(ArgPsi);
            
            if isfield(Arguments,'PsiBC')
                obj.xi0Psi  = Arguments.PsiBC.xi0Psi;
                obj.xi1Psi  = Arguments.PsiBC.xi1Psi;
                obj.xi0PsiTT= Arguments.PsiBC.xi0PsiTT;
            end
            
                                    
%             Arguments.ExtensionParams.Grid          = obj.Grid;
%             Arguments.ExtensionParams.Scatterer     = obj.Scatterer;
%             Arguments.ExtensionParams.Basis         = Arguments.Basis;
%             Arguments.ExtensionParams.CoeffsHandle 	= Arguments.CoeffsHandle;
%             Arguments.ExtensionParams.CoeffsParams	= ArgPsi.CoeffsParams;
%             obj.ExtensionPsi                       = Arguments.Extension(Arguments.ExtensionParams);


        end 
        
		function u = P_Omega(obj,xi_gamma)
			
			rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
			
			rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
			GLW = obj.Gf(rhs(:));
			
			u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));
			u(obj.Scatterer.Np)=xi_gamma(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).' + obj.GF{1}(obj.Scatterer.Np).';
        end

        function P = SPsi(obj,omega)
            
            in = zeros(size(omega));
            in(obj.Scatterer.Mp) = omega(obj.Scatterer.Mp);
            Gw = obj.OpPsi.Solve(in(:));
            
            P=zeros(size(omega));
            
            
            
            u=zeros(size(omega(:)));
            S = obj.Scatterer.TheScatterer();
            u(obj.GridGamma)= obj.xi0Psi(S.th,S.r) + obj.xi1Psi(S.th,S.r).*obj.Scatterer.dr ;
            
            lu=spalloc(size(omega,1),size(omega,2),numel(obj.Scatterer.Mp));
            lu(obj.Scatterer.Mp)=obj.OpPsi.ApplyOp(u,obj.Scatterer.Mp);
            PomegaPsi=u(:) - obj.OpPsi.Solve(lu(:));
            
            
            P(obj.Scatterer.Np) = PomegaPsi(obj.Scatterer.Np) + Gw(obj.Scatterer.Np);
            
        end
      end
      
      methods(Access = protected)
          
          function f = Lu(obj,u,msk)
              if isstruct(u), msk=u.msk; u = u.W; end
              
              if exist('msk','var')
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
          
          function Qj = Qcol(obj,GLW,~)
              Qj = -GLW(obj.GridGamma,:);
          end
          
          function Pc = Pcol(obj,GLW,w)
              
              S = obj.Scatterer.TheScatterer();
              dr = obj.Scatterer.dr;
              dr2=dr.^2;
              
              u=zeros(size(w));
              for i=1:size(w,2)
                u(obj.GridGamma,i)= ...
                            ...%obj.xi0Psi(S.th,S.r) + obj.xi1Psi(S.th,S.r).*dr ...
                            +0.5*(...
                            w(obj.GridGamma,i)...
                            ... %-obj.xi1Psi(S.th,S.r)./S.r ...
                            ... %-obj.xi0PsiTT(S.th,S.r)./(S.r.^2)...
                            ).*dr2;
              end
              %lu=spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
              lu=spalloc(size(w,1),size(w,2),size(w,2)*numel(obj.Scatterer.Mp));
              lu(obj.Scatterer.Mp,:)=obj.OpPsi.ApplyOp(u,obj.Scatterer.Mp);
              glw_psi=obj.OpPsi.Solve(lu);
              
              tmp=zeros(size(w));
              tmp(obj.Scatterer.Mp,:) = w(obj.Scatterer.Mp,:)-GLW(obj.Scatterer.Mp,:);
              %tmp = w-GLW;
              GPsiomega = obj.OpPsi.Solve(tmp);
              Pc = GPsiomega(obj.GridGamma,:) - glw_psi(obj.GridGamma,:);
          end
          
          function GGf = TrGPsiGf(obj,iGf)
              u=zeros(size(iGf));
               S = obj.Scatterer.TheScatterer();
              u(obj.GridGamma)=0.5*iGf(obj.GridGamma).*obj.Scatterer.dr.^2 + obj.xi0Psi(S.th,S.r) + obj.xi1Psi(S.th,S.r).*obj.Scatterer.dr ...
                                              +0.5*(...
                            ...%w(obj.GridGamma,i)...
                            -obj.xi1Psi(S.th,S.r)./S.r ...
                            -obj.xi0PsiTT(S.th,S.r)./(S.r.^2)...
                            ).*obj.Scatterer.dr.^2;
              
              
              lu=spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
              lu(obj.Scatterer.Mp)=obj.OpPsi.ApplyOp(u,obj.Scatterer.Mp);
              glw_psi=obj.OpPsi.Solve(lu(:));

                tGf=zeros(size(iGf));
                tGf(obj.Scatterer.Mp) = iGf(obj.Scatterer.Mp);
               %  tGf=Gf;
                 
              tmp = obj.OpPsi.Solve(tGf);
              GGf = tmp(obj.GridGamma) - glw_psi(obj.GridGamma,:);
          end
 
%           function RhsPsi(obj)not finished
%               obj.ExtensionPsi.Expand();
%               LuPsi = @(u) obj.Op.ApplyOpPsi(u,obj.Scatterer.Mp);
%               
%               tmp=cellfun(@(arg) LuPsi(arg),obj.Extension.W,'UniformOutput',false);
%               
%               obj.rhs = cell(size(tmp));
%               for indx=1:numel(tmp)
%                   [n,m]=size(obj.Extension.W{indx});
%                   NNZ = nnz(obj.Extension.W{indx});
%                   obj.rhs{indx} = spalloc( n,m,NNZ);
%                   obj.rhs{indx}(obj.Scatterer.Mp,:) = tmp{indx};
%               end
%               
%           end
          
          
          function rhs = Bf(obj,F)
              %rhs = F(:);
              rhs = obj.Op.Bf(F);
          end
          
          function res = BF(obj)
              ScattererForSource = obj.Scatterer;
              
              if isequal(obj.SourceHandle , @Tools.Source.NavierStokesDescreteSrc)
                  obj.SourceParams.Np = obj.Scatterer.Np;
                  obj.SourceParams.Mp = obj.Scatterer.Mp;
                  obj.SourceParams.GridGamma = obj.Scatterer.GridGamma;
              end
              
              HS = obj.SourceHandle(ScattererForSource,obj.CoeffsHandle,obj.CoeffsParams,obj.SourceParams);
              
              %res = obj.Bf(HS.Source);
              res = obj.Op.Bf(HS.Source);
              res(obj.Scatterer.Outside())=0;
          end
        end

end

