function NavierStokesTime
    
    LinearSolverType = 3;%0;%
    CollectRhs = 1;
    
    RN = 10; % Reynolds Number
    UsingConvectionTerm=Tools.Enums.Bool.Yes;    
    KindOfConvergance = Tools.Enums.Convergance.Exact;
    
    if 0 %ellipse
        a=1;
        b=0.5;
        
        FocalDistance = sqrt(a^2-b^2);
        Eta0 = acosh(a/FocalDistance);
        
        x1=-1.2;xn=1.2;
        y1=-0.7;yn=0.7;
    else % cirle
        R0 = 1;
        x1=-1.3;xn=1.3;
        y1=-1.3;yn=1.3;
    end
    
    Order=2;
    if Order==2
        Stencil=5;
        Extension = @Tools.Extensions.EBPolarLaplace3OrderExtension;
        %ExtensionPsi = @Tools.Extensions.NavierStokesPsi5rdOrderExtension;
        ExtensionPsi = @Tools.Extensions.NavierStokesPsi4rdOrderExtension;
        
        DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
        %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
    elseif Order==4
        error('doesn''t work')
        Extension = @Tools.Extensions.EBPolarLaplace5OrderExtension;
        
        Stencil=13;
        DiffOp = @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs;
        
    end
    
    ExParams.r0=R0;
    for   r_ = R0*0.9%[0.9,1,1.1]
        ExParams.r = r_;
        
        ExParams.p=99;
        Eparam =@(r) struct('r0',ExParams.r0,'r',r,'p',ExParams.p);
        
        %         BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
        BasisType = Tools.Enums.Basis.Fourier;
        ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
        
        f   =@(phi) Omega(phi,ExParams);
        fn  =@(phi) DrDOmega(phi,ExParams);
        
        Ft = @(t) cossin(t);
        %Ft = @(t) smoothstep(t);

        Basis = BasisType.Helper(f,fn);

                
        ErrOmegaInfPre = 0; ErrOmega2Pre = 0; ErrPsiInfPre=0; ErrPsi2Pre=0; biPre=0; b2Pre=0; pbiPre=0; pb2Pre=0; 
        
        fprintf('Navier Stokes Time:\n Basis - %s, NBss0=%d, NBss1=%d,\n LinearSolverType = %d , \n Order=%d ,Scatterer at r=%-2.2f,Reynold Number=%-2.2f, \n Ft = %s, Convergance:%s,\n UsingConvectionTerm:%s\n', ...
            BasisType.toString(),Basis.NBss0, Basis.NBss1, LinearSolverType,Order,ExParams.r,RN,func2str(Ft),...
            KindOfConvergance.toString(), ...
            UsingConvectionTerm.toString()      );
        firsttime=true;
        
        for n=1:6 %6 %run different grids
            tic
            %build grid
            p=3;%4;
            Nx=2.^(n+p)+1;
            Ny=2.^(n+p)+1;
            
            
            Grid = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
            %ht = abs(Grid.x(1)-Grid.x(2))*0.5;
            ht = Grid.dx;
%tN=1;            
            if firsttime
                tN =round(1.2/ht);
                firsttime=false;
            else
                tN=tN*2;
            end
            k=2*RN/ht;
            
            Fn = @(n) Ft(n*ht);
            
            gn = g0([],ExParams,Grid,k)*Fn(1);
            
            
            SourceParams = struct('Handle',@Tools.Source.NavierStokesSourceRTh,  'SourceParams'      , ExParams,'Src',gn,'Fn',Fn,'r0',ExParams.r0,'ht',ht, 'n',0);
            %gn=zeros(size(gn));
            %gn = g0([],ExParams,Grid,ht)*Fn(0);
            
            xi0Psi  =@(theta,r) Psi(theta,Eparam(r))*Fn(1);
            xi1Psi  =@(theta,r) DrDPsi(theta,Eparam(r))*Fn(1);
            xi0PsiTT=@(theta,r) DPsiDThetaTheta(theta,Eparam(r))*Fn(1);
            xi1PsiTT=@(theta,r) DPsiDrThetaTheta(theta,Eparam(r))*Fn(1);
            xi0PsiTTTT=@(theta,r) DPsiD4Theta(theta,Eparam(r))*Fn(1);
            
            ExtensionParamsPsi.PsiBC = struct('xi0Psi',xi0Psi,'xi1Psi',xi1Psi,'xi0PsiTT',xi0PsiTT,'xi1PsiTT',xi1PsiTT,'xi0PsiTTTT',xi0PsiTTTT);
            
            Setup  = struct( 'Basis', Basis, 'Grid', Grid, ...
                'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                'CoeffsParams'      , struct('a',1,'b',1,'sigma',k), ...
                'ScattererHandle'   , @Tools.Scatterer.PolarScatterer, ...
                'ScattererParams'   , struct('r0',ExParams.r, 'Stencil', Stencil), ...
                'CollectRhs'        , CollectRhs, ...
                'DiffOp'            , DiffOp, ...
                'DiffOpParams'      , struct(   'BC_y1',  0, 'BC_yn',  0,'BC_x1',0 , 'BC_xn',0, 'LinearSolverType', LinearSolverType, 'Order',Order), ...
                'SourceHandle'      , @Tools.Source.NavierStokesDescreteSrc, ...
                'SourceParams'      , SourceParams, ...
                'Extension'         , Extension, ...
                'ExtensionParams'   , [], ...
                'ExtensionPsi'      , ExtensionPsi, ...
                'ExtensionParamsPsi', ExtensionParamsPsi ...
                );
            
            
            Prb =  Solvers.NavierStokesSolver(Setup);
            
            %         CP=Setup.CoeffsParams;
            %         CP.sigma=1;
            %         Src = Tools.Source.NavierStokesSourceRTh(Prb.Scatterer, Setup.CoeffsHandle,CP,ExParams);
            
            % Record = Tools.Common.SimpleRecorder(Grid,tN);
            
            Der = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy);
            
            IsNaN = @(a) any(isnan(full(a(:))));
            IsFinite  = @(a) all(isfinite(full(a(:))));
            
            for tn=1:tN
                %rhs = [-Prb.Qf{1} - Prb.TrGF{1} ; -Prb.TrGGF ];
                %cn = [ Prb.Q{1},Prb.Q{2} ; Prb.P{1},Prb.P{2}]\rhs;
                rhs = [-Prb.Qf{1}-Prb.TrGF{1} ; -Prb.Qpsi{1} - Prb.QPsif - Prb.TrGGF - Prb.TrGpsiGExf ];
                cn = [ Prb.Q{1},Prb.Q{2} ; Prb.QpsiOmega{1} + Prb.GPO{1}, Prb.QpsiOmega{2} + Prb.GPO{2}]\rhs;
                
                Oxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                Oxi(Prb.GridGamma) = [Prb.W{1}(Prb.GridGamma,:),Prb.W{2}(Prb.GridGamma,:)]*cn  + Prb.Wf{1}(Prb.GridGamma) + Prb.WPsif{1}(Prb.GridGamma);
                
                uOmega = Prb.P_Omega(Oxi);
                assert(~IsNaN(uOmega))
                
                %if UsingConvectionTerm
                    
                    Pxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                    Pxi(Prb.GridGamma) = [Prb.WpsiOmega{1}(Prb.GridGamma,:),Prb.WpsiOmega{2}(Prb.GridGamma,:)]*cn + Prb.Wpsi{1}(Prb.GridGamma,:);
                    
                    uPsi = Prb.P_OmegaPsi(uOmega,Pxi);
                    assert(~IsNaN(uPsi))
                    
                 if UsingConvectionTerm   
                    [Omega_x,Omega_y]       = Der.CartesianDerivatives(uOmega);
                    [Psi_x,Psi_y]       = Der.CartesianDerivatives(uPsi);
                    
                    D = zeros(Nx,Ny); 
                    tmp = Psi_y.*Omega_x - Psi_x.*Omega_y;
                    D(Prb.Scatterer.Mp) = tmp(Prb.Scatterer.Mp);
                else
                    D = 0;
                end
                
               assert(IsFinite(D) )
                
                % Record.Store(tn*ht,u);
                
                % Prepare for next time step
                
                %ftilde = -Src.Source()*(Fn(tn)+Fn(tn+1));
                ftilde = Ftilde(ExParams,Grid,Ft,ht*tn,RN,UsingConvectionTerm) + Ftilde(ExParams,Grid,Ft,ht*(tn+1),RN,UsingConvectionTerm);
                
                
                
                gn = 2*RN*D -(gn + 4*RN*uOmega/ht +  RN*ftilde);
                SourceParams.Src = gn ;
                SourceParams.n = tn;
                SourceParams.ht = ht;
                UpdateParams.SourceParams = SourceParams;
                
                UpdateParams.PsiBC.xi0Psi       = @(theta,r) Psi(theta,Eparam(r))*Fn(tn+1);
                UpdateParams.PsiBC.xi1Psi       = @(theta,r) DrDPsi(theta,Eparam(r))*Fn(tn+1);
                UpdateParams.PsiBC.xi0PsiTT     = @(theta,r) DPsiDThetaTheta(theta,Eparam(r))*Fn(tn+1);
                UpdateParams.PsiBC.xi1PsiTT     = @(theta,r) DPsiDrThetaTheta(theta,Eparam(r))*Fn(tn+1);
                UpdateParams.PsiBC.xi0PsiTTTT   = @(theta,r) DPsiD4Theta(theta,Eparam(r))*Fn(tn+1);
                
                Prb.Update( UpdateParams );
            end
            
            %save('Movie.mat','Record');
            
            if 0% ~UsingConvectionTerm

                Pxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                Pxi(Prb.GridGamma) = [Prb.WpsiOmega{1}(Prb.GridGamma,:),Prb.WpsiOmega{2}(Prb.GridGamma,:)]*cn + Prb.Wpsi{1}(Prb.GridGamma,:);

                uPsi = Prb.P_OmegaPsi(uOmega,Pxi);
            end
            t1=toc;
            
            %------------------------------------------------------------------
            % Comparison
            %------------------------------------------------------------------
            
            
            
            if strcmpi(KindOfConvergance,'Grid')
                if n > 1
                    uOmega1= spalloc(Nx,Ny,nnz(uOmega0));
                    uOmega1(Prb.Np) = uOmega0(Prb.Np);
                    
                    tmp = uOmega(1:2:end,1:2:end)-uOmega1(1:2:end,1:2:end);
                    ErrOmegaInf = norm(tmp(:),inf);
                    
                    
                    uPsi1= spalloc(Nx,Ny,nnz(uPsi0));
                    uPsi1(Prb.Np) = uPsi0(Prb.Np);
                    
                    tmp = uPsi(1:2:end,1:2:end)-uPsi1(1:2:end,1:2:end);
                    ErrPsiInf = norm(tmp(:),inf);
                    
                    fprintf('N=%-10dx%-10d\t ErrTot=%d\t rate=%-5.2f\t ErrTot=%d\t rate=%-5.2f\t time=%d\n', Nx,Ny, ...
                        ErrOmegaInf,log2(ErrOmegaInfPre/ErrOmegaInf),ErrPsiInf,log2(ErrPsiInfPre/ErrPsiInf),t1);
                    
                    ErrOmegaInfPre = ErrOmegaInf;
                    ErrPsiInfPre = ErrPsiInf;
                                        
                end
                
                uOmega0=spalloc(Nx*2-1,Ny*2-1,nnz(uOmega));
                uOmega0(1:2:end,1:2:end)=uOmega;
                
                uPsi0=spalloc(Nx*2-1,Ny*2-1,nnz(uPsi));
                uPsi0(1:2:end,1:2:end)=uPsi;

                
            elseif strcmpi(KindOfConvergance,'Exact')
                ExParams2 = ExParams;
                ExParams2.r = Prb.Scatterer.r;
                OxiEx = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                OxiEx(Prb.GridGamma) = Omega(Prb.Scatterer.th,ExParams2)*Fn(tn);
                [bi,b2] = cmpr(OxiEx(Prb.GridGamma), Oxi(Prb.GridGamma));
                
                PxiEx = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                PxiEx(Prb.GridGamma) = Psi(Prb.Scatterer.th,ExParams2)*Fn(tn);
                [pbi,pb2] = cmpr(PxiEx(Prb.GridGamma), Pxi(Prb.GridGamma));
                
                
                Ex = zeros(size(Grid.R));
                
                %Ex(IntPrb.Scatterer.Np) = Ex(FocalDistance,Prb.Scatterer.R(Prb.Scatterer.Np),Prb.Scatterer.Th(Prb.Scatterer.Np),R0);
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.Scatterer.Np);
                Ex(Prb.Scatterer.Np) = Omega(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3)*Fn(tn);
                
                %MakeExactMovie(Grid,tN,Prb.Scatterer.Np,Prb.Scatterer.Th,ExParams3,Fn,ht);
                
                ExP = zeros(size(Grid.R));
                ExP(Prb.Scatterer.Np) = Psi(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3)*Fn(tn);
                
                msk = Prb.Scatterer.Np;
                [ErrOmegaInf,ErrOmega2] = cmpr(Ex(msk),uOmega(msk));
                [ErrPsiInf,ErrPsi2]     = cmpr(ExP(msk),uPsi(msk));
                
                str = sprintf('%-8.0f %-8.4f %-8.1f %-8.4f $%-6d\\\\times %-7d$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-6.2f$ \\\\\\\\ \n',...
                    k,ht,tN,tn*ht,...
                    Nx,Ny,ErrOmegaInf, log2(ErrOmegaInfPre/ErrOmegaInf), ...
                    ErrPsiInf,  log2(ErrPsiInfPre/ErrPsiInf)  , ...
                    bi,  log2(biPre/bi) , ... 
                    pbi,  log2(pbiPre/pbi) , ...
                    t1         );
                
                ErrOmegaInfPre = ErrOmegaInf; ErrOmega2Pre = ErrOmega2;
                biPre=bi; b2Pre=b2;
                pbiPre=pbi; pb2Pre=pb2;                
                ErrPsiInfPre = ErrPsiInf; ErrPsi2Pre = ErrPsi2;
                
                fprintf(str);
                %fprintf(fileID,str);
                 
            end            
        end
    end
end


function [Linf,L2] = cmpr(ex,u)%,GG)
    
    tmp = ex - u;
    
    %     if exist('GG','var')
    %         tmp(GG) = 0;
    %         u(GG)=0;
    %     end
    
    Linf = norm(tmp(:),inf);%/norm(u(:),inf);
    L2   = norm(tmp(:),2);%/norm(u(:),2);
end

function P = Psi(theta,Params)
    %         x = r .* cos(theta);
    %         y = r .* sin(theta);
    assert(Params.p==99);
    
    P = (Params.r.^2).*cos(theta).*sin(theta).*(Params.r.^2-Params.r0.^2).^2;%Params.p;
end

function DP = DrDPsi(theta,Params)
    %         x = r .* cos(theta);
    %         y = r .* sin(theta);
    
    DP = Params.r.*(Params.r.^2-Params.r0.^2).*(3*Params.r.^2-Params.r0.^2).*sin(2*theta);
end

function DP = DPsiDThetaTheta(theta,Params)
    %         x = r .* cos(theta);
    %         y = r .* sin(theta);
    
    DP = -4*Psi(theta,Params);
    
end

function DP = DPsiD4Theta(theta,Params)
    %         x = r .* cos(theta);
    %         y = r .* sin(theta);
    
    DP = 16*Psi(theta,Params);
    
end

function DP = DPsiDrThetaTheta(theta,Params)
    %         x = r .* cos(theta);
    %         y = r .* sin(theta);
    
    DP = -4*DrDPsi(theta,Params);
    
end

function O = Omega(theta,Params)
    
    p=Params.p;
    r0=Params.r0;
    if numel(Params.r)==1
        r=ones(size(theta))*Params.r;
    else
        r = Params.r;
    end
    
    O =  4*(r.^2).*(4*r.^2-3*r0^2).*sin(2*theta);
end

function DO = DrDOmega(theta,Params)
    assert(Params.p>=2);
    
    p=Params.p;
    r0=Params.r0;
    
    if numel(Params.r)==1
        r=ones(size(theta))*Params.r;
    else
        r = Params.r;
    end
    
    switch p
        case 0
            DO=0;
            warning('with p=0 you get Omega=0, so Psi does not satisfy the BC');
        case 1
            DO=0;
            
        case 2
            DO = 32*r.^2;
        case 99 %special case (r, theta)
            DO=8*r.*(8*r.^2-3*r0^2).*sin(2*theta);
        otherwise
            DO = 8*(p-1)*p.*r.*((r.^2 - r0^2).^(p-3)).*(p.*r.^2 - 2*r0^2);
    end
    
end


function MakeExactMovie(Grid,tN, Np,Th,Params,Fn,ht)
    
    Record = Tools.Common.SimpleRecorder(Grid,tN);
    
    Ex = zeros(size(Grid.R));
    
    for tn=1:tN
        Ex(Np) = Omega(Th(Np),Params)*Fn(tn);
        
        Record.Store(tn*ht,Ex);
    end
    
    save('ExOmegaCos.mat','Record');
    
end

function [y,dy] = cossin(t)
    y = cos(t);
    dy = -sin(t);
end

function g=g0(~,Params,Grid,k)
    
    Params.r = Grid.R;
    g = 384*Grid.X.*Grid.Y ... - 8*Grid.X.*Grid.Y.*(4*(Grid.X.*Grid.X+Grid.Y.*Grid.Y) - 3*Params.r0*Params.r0)*(2/h);
        -k*Omega(Grid.Theta,Params);
    
    % -8*Grid.X.*Grid.Y.*(48+(4*r.^2-3*r0^2));
    
    %g=zeros(size(g));
end

function g=Ftilde(Params,Grid,Func,t,RN,UseConvTerm)
    
    x=Grid.X;
    y=Grid.Y;
    xy = x.*y;
    r = Grid.R;
    
    Params.r = r;
    [ft,dft] = Func(t);
    
    g = Omega(Grid.Theta,Params)*dft - 384*xy*ft/RN;
    
    if UseConvTerm
        D= -64*xy.*(x.^2-y.^2).*((r.^2 - 3*(Params.r0^2)/4).^2 - (Params.r0^4)/16);
        g=g+D.*ft.*ft;
    end
end