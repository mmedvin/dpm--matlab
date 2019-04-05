classdef NavierStokesSolver < Solvers.SuperNoNHomoNavierStokesSolver
    %NAVIERSTOCKSHOMOSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
%     properties(Access = protected, AbortSet = true)
%         Op;
%     end
%     properties(Access = public)%???
%         
%     end
    
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
            obj.xi0PsiTT= Params.PsiBC.xi0PsiTT
            obj.xi0PsiTTTT= Arguments.PsiBC.xi0PsiTTTT;
            obj.xi1PsiTT= Params.PsiBC.xi1PsiTT;
            
            obj.CreateWf();
            obj.calc_QnWf();          
        end
        
        function obj = NavierStokesSolver(Arguments)
            %SUPERHOMONAVIERSTOCKSSOLVER Construct an instance of this class
            %   Detailed explanation goes here
            
            obj = obj@Solvers.SuperNoNHomoNavierStokesSolver(Arguments);
            
            

            
            
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
        
        function u = P_OmegaPsi(obj,omega,xi_omegaPsi)
            
             rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
             
             rhs(obj.Scatterer.Mp)= obj.Lu(xi_omegaPsi(:),obj.Scatterer.Mp,obj.OpPsi);
             GLW = obj.Gf(rhs(:),obj.OpPsi);

%              lu=spalloc(size(omega,1),size(omega,2),numel(obj.Scatterer.Mp));
%              lu(obj.Scatterer.Mp)=obj.OpPsi.ApplyOp(xi_omegaPsi(:),obj.Scatterer.Mp);
%              PomegaPsi=xi_omegaPsi(:) - obj.OpPsi.Solve(lu(:));


            in = zeros(size(omega));
            in(obj.Scatterer.Mp) = omega(obj.Scatterer.Mp);
            Gw = obj.OpPsi.Solve(in(:));

            u = spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Np));                        
            u(obj.Scatterer.Np) =  xi_omegaPsi(obj.Scatterer.Np) - GLW(obj.Scatterer.Np).' ...
            ... PomegaPsi(obj.Scatterer.Np)
            + Gw(obj.Scatterer.Np).';
            

            
        end        
        
        function P = SPsiDeleteMe(obj,omega,Pxi)%O,Or,Ott)
            
            in = zeros(size(omega));
            in(obj.Scatterer.Mp) = omega(obj.Scatterer.Mp);
            Gw = obj.OpPsi.Solve(in(:));
            
          %  P=zeros(size(omega));
            
            
            
            %u=zeros(size(omega(:)));
            %S = obj.Scatterer.TheScatterer();
%             dr = obj.Scatterer.dr;
%             dr2=dr.^2;
%             dr3=dr2.*dr;
%             dr4=dr3.*dr;
%             
%             r = S.r;
%             r2=r.*r;
%             r3=r2.*r;
%             r4=r3.*r;
            
%             xi0P = obj.xi0Psi(S.th,S.r);
%             xi1P = obj.xi1Psi(S.th,S.r);
%             xi0Ptt = obj.xi0PsiTT(S.th,S.r);
%             xi0Ptttt = obj.xi0PsiTTTT(S.th,S.r);
%             xi1Ptt = obj.xi1PsiTT(S.th,S.r);
            
%             Orr = obj.Wf{1}(obj.GridGamma) - Or./r - Ott./r2 + obj.CoeffsParams.sigma*O;
%             
%             Prr = O - xi1P./r - xi0Ptt./r2;
%             Prrtt = Ott - xi1Ptt./r - xi0Ptttt./r2;
%             
%             P3r = Or + xi1P./r2 - Prr./r + 2*xi0Ptt./r3 - xi1Ptt./r2;
%             P4r = Orr - 2*xi1P./r3 + 2*Prr./r2 - P3r./r - 6*xi0Ptt./r4 + 4*xi1Ptt./r3 - Prrtt./r2;
%             
%             u(obj.GridGamma)= xi0P + xi1P.*dr + Prr.*(dr2/2)+ P3r.*(dr3/6) + P4r.*(dr4/24);
            
            %u(obj.GridGamma) = obj.ExtendPsi(S.r,O,Or,Ott,xi0P,xi1P,xi0Ptt,xi1Ptt,xi0Ptttt,obj.Wf{1}(obj.GridGamma));
            
            
            %u(obj.GridGamma)= obj.xi0Psi(obj.Scatterer.th,obj.Scatterer.r);%debug
            
%             lu=spalloc(size(omega,1),size(omega,2),numel(obj.Scatterer.Mp));
%             lu(obj.Scatterer.Mp)=obj.OpPsi.ApplyOp(u,obj.Scatterer.Mp);
%             PomegaPsi=u(:) - obj.OpPsi.Solve(lu(:));
%             
            rhs=zeros(obj.Grid.Nx,obj.Grid.Ny);
            
            rhs(obj.Scatterer.Mp)= obj.Lu(Pxi(:),obj.Scatterer.Mp,obj.OpPsi);
            GLW = obj.Gf(rhs(:),obj.OpPsi);
            PomegaPsi = Pxi(obj.Scatterer.Np) - GLW(obj.Scatterer.Np);
            
            P(obj.Scatterer.Np) =  PomegaPsi +  Gw(obj.Scatterer.Np);
            
        end
    end
    
    methods(Access = protected)
        
        function f = Lu(obj,u,msk,OP)
            if isstruct(u), msk=u.msk; u = u.W; end
            if ~exist('OP','var'), OP = obj.Op; end
            
            if exist('msk','var')
                %f = obj.Op.A(msk,:)*u;
                f = OP.ApplyOp(u,msk);
            else
                %f = obj.Op.A*u;
                f = OP.ApplyOp(u);
            end
        end
        
        function u = Gf(obj,f,OP)
            if ~exist('OP','var'), OP = obj.Op; end
            
            %u = obj.Op.A\f;
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
        
%         function Pc = Pcol(obj,GLW,w)
%             
% %             S = obj.Scatterer.TheScatterer();
% %             dr = obj.Scatterer.dr;
% %             dr2=dr.^2;
% %             dr3=dr2.*dr;
% %                         
% %             u=zeros(size(w)); %Tr u = xiPsi (the part with omega)
% %             
% %             %u(obj.GridGamma,1) = (dr2/2 - dr3/6./S.r).*w(obj.GridGamma,1);
% %             %u(obj.GridGamma,2) = (dr3/6).* w(obj.GridGamma,2);
% %             
% %             %here I use overuse the fact about known exact Omegam works only with specific one
% %             Ott = -4*w(obj.GridGamma,1);
% %             u(obj.GridGamma,1) = obj.ExtendPsi(S.r,w(obj.GridGamma,1),0                 ,Ott,0,0,0,0,0,0);
% %             u(obj.GridGamma,2) = obj.ExtendPsi(S.r,0                 ,w(obj.GridGamma,2),0  ,0,0,0,0,0,0);
% %             
% %             
% %             % computing potential
% %             Lu=spalloc(size(w,1),size(w,2),size(w,2)*numel(obj.Scatterer.Mp));
% %             Lu(obj.Scatterer.Mp,:)=obj.OpPsi.ApplyOp(u,obj.Scatterer.Mp);
% %             glw_psi= - obj.OpPsi.Solve(Lu);% (P-I)_\gamma^{(\psi)}, i.e. the BEP
%             
%             PO=zeros(size(w)); % P_\Omega^{(\omega)} \xi_\gamma^{(\omega)}
%             PO(obj.Scatterer.Mp,:) = w(obj.Scatterer.Mp,:)-GLW(obj.Scatterer.Mp,:);
%             
%             %GPsiomega is actually G^{(\psi)} PO, however the actual G^{(\psi)} \omega  is linear combination of  G^{(\psi)} PO
%             GPsiomega = obj.OpPsi.Solve(PO);
%             
%             % (P-I)_\gamma^{(\psi)} +  G^{(\psi)} PO
%             Pc =  GPsiomega(obj.GridGamma,:); %+glw_psi(obj.GridGamma,:) ;
%         end
%         
        function GGf = TrGPsiGf(obj,iGf)
            
            S = obj.Scatterer.TheScatterer();
            dr = obj.Scatterer.dr;
            dr2=dr.^2;
            dr3=dr2.*dr;
 
 
            u=zeros(size(iGf));
            
           % u(obj.GridGamma)= obj.xi0Psi(S.th,S.r) ...
           %                 +(dr - (dr2/2)./S.r + 2*(dr3/6)./(S.r.^2)).*obj.xi1Psi(S.th,S.r) ...
           %                 -((dr2/2)./(S.r.^2) - 3*(dr3/6)./(S.r.^3)).*obj.xi0PsiTT(S.th,S.r) ...
           %                 -((dr3/6)./(S.r.^2)).*obj.xi1PsiTT(S.th,S.r);
           

%            Source = SourceHandle(obj.Scatterer.TheScatterer,obj.CoeffsHandle,obj.CoeffsParams,SourceParams);
%            NoXiPsi = Tools.Misc.XiPsi();
%            u(obj.GridGamma) = obj.ExtendPsi(S.r,NoXiPsi,NoXiPsi,obj.xiPsi,Source);
%            
%             
%             lu=spalloc(obj.Grid.Nx,obj.Grid.Ny,numel(obj.Scatterer.Mp));
%             lu(obj.Scatterer.Mp)=obj.OpPsi.ApplyOp(u,obj.Scatterer.Mp);
%             glw_psi=-obj.OpPsi.Solve(lu(:));
            
            tGf=zeros(size(iGf));
            tGf(obj.Scatterer.Mp) = iGf(obj.Scatterer.Mp);
            %  tGf=Gf;
            
            tmp = obj.OpPsi.Solve(tGf);
            GGf = tmp(obj.GridGamma);% + glw_psi(obj.GridGamma,:);
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

