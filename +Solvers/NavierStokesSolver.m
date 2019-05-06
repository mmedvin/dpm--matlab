classdef NavierStokesSolver < Solvers.SuperNoNHomoNavierStokesSolver
    %NAVIERSTOCKSHOMOSOLVER 
    
    methods
        
        function Update(obj,Params)
            obj.UpdatePsi(Params.PsiBC);
            obj.UpdateSource(Params);
        end
        
        function obj = NavierStokesSolver(Arguments)
            obj = obj@Solvers.SuperNoNHomoNavierStokesSolver(Arguments);         
        end
        
        function u = P_Omega(obj,xi_gamma)
            
            rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
            
            rhs(obj.Scatterer.Mp)= obj.Lu(xi_gamma(:),obj.Scatterer.Mp);
            GLW = obj.Gf(rhs(:));
            
            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));
            u(obj.Scatterer.Np)=xi_gamma(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).' + obj.GF{1}(obj.Scatterer.Np).';
        end
        
        function u = P_OmegaPsi(obj,omega,xi_omegaPsi)
            
             rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
             
             rhs(obj.Scatterer.Mp)= obj.Lu(xi_omegaPsi(:),obj.Scatterer.Mp,obj.OpPsi);
             GLW = obj.Gf(rhs(:),obj.OpPsi);

            in = zeros(size(omega));
            in(obj.Scatterer.Mp) = omega(obj.Scatterer.Mp);
            Gw = obj.OpPsi.Solve(in(:));

            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));                        
            u(obj.Scatterer.Np) =  xi_omegaPsi(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).' + Gw(obj.Scatterer.Np).';
                        
        end                
    end
    
    methods(Access = protected)
        
        function f = Lu(obj,u,msk,OP)
            if isstruct(u), msk=u.msk; u = u.W; end
            if ~exist('OP','var'), OP = obj.Op; end
            
            if exist('msk','var')
                f = OP.ApplyOp(u,msk);
            else
                f = OP.ApplyOp(u);
            end
        end
        
        function u = Gf(obj,f,OP)
            if ~exist('OP','var'), OP = obj.Op; end
            
            if numel(f) == obj.Grid.Nx*obj.Grid.Ny
                u = OP.Solve(f(:));
            else
                u = OP.Solve(f);
            end
        end
        
        function Qj = Qcol(obj,GLW,~)
            Qj = -GLW(obj.GridGamma,:);
        end
        
        function TGP = TrGpsiPOmega(obj,GLW,w)
        
                   
            PO=zeros(size(w)); % P_\Omega^{(\omega)} \xi_\gamma^{(\omega)}
            PO(obj.Scatterer.Mp,:) = w(obj.Scatterer.Mp,:)-GLW(obj.Scatterer.Mp,:);
            
            %GPsiomega is actually G^{(\psi)} PO, however the actual G^{(\psi)} \omega  is linear combination of  G^{(\psi)} PO
            GPsiomega = obj.Gf(PO,obj.OpPsi); %obj.OpPsi.Solve(PO);%
            
            % (P-I)_\gamma^{(\psi)} +  G^{(\psi)} PO
            TGP =  GPsiomega(obj.GridGamma,:); %+glw_psi(obj.GridGamma,:) ;

        end
                
        function GGf = TrGPsiGf(obj,iGf)
            
            S = obj.Scatterer.TheScatterer();
            dr = obj.Scatterer.dr;
            dr2=dr.^2;
            dr3=dr2.*dr;
 
 
            u=zeros(size(iGf));
                        
            tGf=zeros(size(iGf));
            tGf(obj.Scatterer.Mp) = iGf(obj.Scatterer.Mp);
            %  tGf=Gf;
            
            tmp = obj.OpPsi.Solve(tGf);
            GGf = tmp(obj.GridGamma);% + glw_psi(obj.GridGamma,:);
        end
               
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

