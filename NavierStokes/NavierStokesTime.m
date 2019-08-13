function NavierStokesTime
    
    ExactChoice = NavierStokesExact.Exact1Time;
    Exact = ExactChoice.Helper();
    
    RN = 10; % Reynolds Number
    UsingConvectionTerm     = Tools.Enums.Bool.No;
    UsingNumericalLaplacian = Tools.Enums.Bool.No;
    
    Order=2;
    [SetupBase, BasisType, ExParams,GridParams,KindOfConvergance] = Params(Order);
    
    Ft =  @cossin;
    %Ft = @smoothstep;
    
    RecordMovie = false;
    nToRecord = 5;
    
    for   r_ = ExParams.r0*0.9%[0.9,1,1.1]
        ExParams.r = r_;
        Eparam =@(r) struct('r0',ExParams.r0,'r',r);
        
        f   =@(phi) Exact.Omega(phi,ExParams);
        fn  =@(phi) Exact.DrDOmega(phi,ExParams);
        Basis = BasisType.Helper(f,fn);
        
        ErrOmegaInfPre = 0; ErrOmega2Pre = 0; ErrPsiInfPre=0; ErrPsi2Pre=0; biPre=0; b2Pre=0; pbiPre=0; pb2Pre=0;
        
        fprintf('Navier Stokes Time:\n Basis - %s, NBss0=%d, NBss1=%d, LinearSolverType = %d , \n Order=%d ,Scatterer at r=%-2.2f,Reynold Number=%-2.2f, \n', ...
            BasisType.toString(),Basis.NBss0, Basis.NBss1,SetupBase.DiffOpParams.LinearSolverType,Order,ExParams.r,RN);
        
        fprintf('%s, Convergance:%s,\n UsingConvectionTerm:%s\n', ...
            ExactChoice.toString(Ft), KindOfConvergance.toString(), UsingConvectionTerm.toString()      );
 
        
        firsttime=true;        
        for n=1:6 %6 %run different grids
            tic
            %build grid
            p=3;%4;
            [Nx,Ny] = deal(2.^(n+p)+1,2.^(n+p)+1);
             
            UpdateParams.SourceParams.Src = zeros(Nx,Ny);
            Grid = Tools.Grid.CartesianGrid(GridParams.x1,GridParams.xn,Nx,GridParams.y1,GridParams.yn,Ny);
            ht = Grid.dx;
            if firsttime
                tN =round(1.1/ht);
                firsttime=false;
            else
                tN=tN*2;
            end
            k=2*RN/ht;
            
            Fn = @(n) Ft(n*ht);
            
            %TStpSrc = @(theta,r,gn,Omega,n,D) TimeStepSource(theta,Eparam(r),gn,Omega,Ft,ht,n,RN,UsingConvectionTerm,Exact,D);
                        
                                                
            Setup  = SetupBase; Setup.Basis = Basis; Setup.Grid = Grid;
            Setup.CoeffsParams.sigma = k;
            Setup.ScattererParams.r0 = ExParams.r;
            
            %gn = Exact.L_np1(Grid.Theta,Eparam(Grid.R), k,Fn,1); %g1 
            gn = g1(ExParams,k,Fn,Grid,Setup.ScattererParams,Exact);

            
            t1=1;
            Setup.SourceParams = struct('Src',gn,'Fn',Fn,'r0',ExParams.r0,'ht',ht, 'tn',t1,'RN',RN);%,'TStpSrc',TStpSrc,'Exact',Exact,'k',k);
            
            Setup.ExtensionParamsPsi.PsiBC = PsiBC(Eparam,Exact,Fn,t1);
                        
            Prb =  Solvers.NavierStokesSolver(Setup);
            
            if RecordMovie  %initialize
                RecorderOmega = Tools.Common.SimpleRecorder(Grid,tN);
                RecorderPsi = Tools.Common.SimpleRecorder(Grid,tN);
            end
            if UsingConvectionTerm || UsingNumericalLaplacian, Der = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy); end
            
            %for use in each time step
            Zeros =zeros(Grid.Nx,Grid.Ny);
            gnp1 = Zeros;
            uO = Zeros;
            Msk = Prb.Scatterer.Mp;
            theta = Zeros;
            theta(Msk) = Grid.Theta(Msk);
            rr=Zeros;
            rr(Msk) = Grid.R(Msk);
            
            for tn=t1:tN
                rhs = [-Prb.Qf{1}-Prb.TrGF{1} ; -Prb.Qpsi{1} - Prb.QPsif - Prb.TrGGF - Prb.TrGpsiGExf ];
                cn = [ Prb.Q{1},Prb.Q{2} ; Prb.QpsiOmega{1} + Prb.GPO{1}, Prb.QpsiOmega{2} + Prb.GPO{2}]\rhs;
                
                Oxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                Oxi(Prb.GridGamma) = [Prb.W{1}(Prb.GridGamma,:),Prb.W{2}(Prb.GridGamma,:)]*cn  + Prb.Wf{1}(Prb.GridGamma) + Prb.WPsif{1}(Prb.GridGamma);
                
                %Oxi_=Oxi;
                %Oxi_(Prb.Scatterer.GridGamma)  = Exact.Omega(Grid.Theta(Prb.Scatterer.GridGamma),Eparam(Grid.R(Prb.Scatterer.GridGamma)))*Fn(tn);
                uOmega = Prb.P_Omega(Oxi);
                
                if RecordMovie && n==nToRecord, RecorderOmega.Store((tn-1)*ht,uOmega);  end
                
                if tn==tN || UsingConvectionTerm
                    
                    Pxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                    Pxi(Prb.GridGamma) = [Prb.WpsiOmega{1}(Prb.GridGamma,:),Prb.WpsiOmega{2}(Prb.GridGamma,:)]*cn + Prb.Wpsi{1}(Prb.GridGamma,:);
                    
                    uPsi = Prb.P_OmegaPsi(uOmega,Pxi);
                end
                
                if UsingConvectionTerm
                    RecorderPsi.Store(tn*ht,uOmega);
                    
                    [Omega_x,Omega_y]   = Der.CartesianDerivatives(uOmega);
                    [Psi_x,Psi_y]       = Der.CartesianDerivatives(uPsi);
                    D = Psi_y.*Omega_x - Psi_x.*Omega_y;
                else
                    D = 0;
                end
                
                if UsingNumericalLaplacian
                    LapO = Der.CartesianLaplacian(uOmega);
                                        
                    Src = @(n) Exact.NSTimeSource(Grid.Theta,Eparam(Grid.R),RN,Ft,tn*ht,UsingConvectionTerm);

                    gn = Exact.DL_n(LapO - k*uOmega,uOmega,Grid.Theta,Eparam(Grid.R), RN,ht,Ft,tn,UsingConvectionTerm);

                    %gn = 2*RN*D -(LapO + k*uOmega +  RN*(Src(tn+1) + Src(tn)) );
                    
                else
                    if 0     
                        uO(Msk) = uOmega(Msk);  
                         
                        gnp1 = Exact.DL_n(gn,uO,theta,Eparam(rr), RN,ht,Ft,tn,UsingConvectionTerm);
                        gn = gnp1;
                    else
                        uO=zeros(size(uOmega));
                        Msk = Prb.Scatterer.Mp;
                        uO(Msk)  = Exact.Omega(Grid.Theta(Msk),Eparam(Grid.R(Msk)))*Fn(tn);
                        gnpre = gn;
                        %norm(uOmega(:)-uO(:),inf)
                        
                        gn=zeros(size(uOmega)); 
                        Msk = Prb.Scatterer.Mp;
                        gn(Msk) = Exact.DL_n(gnpre(Msk),uO(Msk)    ,Grid.Theta(Msk),Eparam(Grid.R(Msk)), RN,ht,Ft,tn,UsingConvectionTerm);
                    end

                    %gn2=gn;
                    %gn2(Msk)= Exact.DL_n(gn(Msk),uOmega(Msk),Grid.Theta(Msk),Eparam(Grid.R(Msk)), RN,ht,Ft,tn,UsingConvectionTerm);
                    %norm(gn(Msk)-gn2(Msk),inf)
                    
                    %gn = TimeStepSource(Grid.Theta,Eparam(Grid.R),pregn,uOmega,Ft,ht,n,RN,UsingConvectionTerm,Exact,D);
                    
                    %gn = TStpSrc(Grid.Theta,Grid.R,pregn,uOmega,tn,D);
                    
                   % gn = Exact.L_np1(Grid.Theta,Eparam(Grid.R), k,Fn,tn+1);
                    %Res = Exact.Test(Grid.Theta,Eparam(Grid.R),RN,ht,Ft,tn);
                    
                end
            
                UpdateParams.SourceParams = Setup.SourceParams;
                UpdateParams.SourceParams.Src = gn;%(Prb.Scatterer.Mp) = gn(Prb.Scatterer.Mp) ;
                UpdateParams.SourceParams.tn   = tn+1;% n+1
                UpdateParams.SourceParams.ht  = ht;
                
                UpdateParams.PsiBC = PsiBC(Eparam,Exact,Fn,tn+1);
                                
                Prb.Update( UpdateParams );
                
            end
            
            if RecordMovie
                save('MovieOmega.mat','RecordOmega');
                if UsingConvectionTerm
                    save('MoviePsi.mat','RecordPsi');
                end
            end
            
            %if n==nToRecord, MakeExactMovie(Grid,tN, Prb.Np,ExParams,Ft,ht, Exact);  end
            
            
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
                OxiEx(Prb.GridGamma) = Exact.Omega(Prb.Scatterer.th,ExParams2)*Fn(tN);
                [bi,b2] = cmpr(OxiEx(Prb.GridGamma), Oxi(Prb.GridGamma));
                
                PxiEx = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                PxiEx(Prb.GridGamma) = Exact.Psi(Prb.Scatterer.th,ExParams2)*Fn(tN);
                [pbi,pb2] = cmpr(PxiEx(Prb.GridGamma), Pxi(Prb.GridGamma));
                
                
                Ex = zeros(size(Grid.R));
                
                %Ex(IntPrb.Scatterer.Np) = Ex(FocalDistance,Prb.Scatterer.R(Prb.Scatterer.Np),Prb.Scatterer.Th(Prb.Scatterer.Np),R0);
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.Scatterer.Np);
                Ex(Prb.Scatterer.Np) = Exact.Omega(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3)*Fn(tN);
                
                %MakeExactMovie(Grid,tN,Prb.Scatterer.Np,Prb.Scatterer.Th,ExParams3,Fn,ht);
                
                ExP = zeros(size(Grid.R));
                ExP(Prb.Scatterer.Np) = Exact.Psi(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3)*Fn(tN);
                
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

function [SetupBase, BasisType, ExParams,GridParams,KindOfConvergance] = Params(Order)
    if nargin==0, Order=2;end
    
    LinearSolverType = 3;%0;%
    CollectRhs = 1;
    
    KindOfConvergance = Tools.Enums.Convergance.Exact;
    BasisType = Tools.Enums.Basis.Fourier;
    ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
    
    
    if Order==2
        Stencil=5;
        Extension = @Tools.Extensions.EBPolarSimpleLaplace3OrderExtension;
        %Extension = @Tools.Extensions.EBPolarLaplace3OrderExtension;
        %Extension = @Tools.Extensions.EBPolarLaplace5OrderExtension;
        ExtensionPsi = @Tools.Extensions.NavierStokesPsi5rdOrderExtension;
        %ExtensionPsi = @Tools.Extensions.NavierStokesPsi4rdOrderExtension;
        
        DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
        %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
        
        DiffOpParams = struct(   'BC_y1',  0, 'BC_yn',  0,'BC_x1',0 , 'BC_xn',0, 'LinearSolverType', LinearSolverType, 'Order',Order);
    elseif Order==4
        error('doesn''t work')
        Extension = @Tools.Extensions.EBPolarLaplace5OrderExtension;
        
        Stencil=13;
        DiffOp = @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs;
        
    end
    
    if 1         % cirle
        ExParams.r0 = 1;
        GridParams.x1=-1.3;
        GridParams.xn=1.3;
        GridParams.y1=-1.3;
        GridParams.yn=1.3;
    else  %ellipse
        a=R0;
        b=a/2;
        
        FocalDistance = sqrt(a^2-b^2);
        Eta0 = acosh(a/FocalDistance);
        
        x1=-1.2;xn=1.2;
        y1=-0.7;yn=0.7;
    end
    
    
    
    SetupBase  = struct( ...
        ... % 'Basis'            , Basis, ...
        ... % 'Grid'             , Grid, ...
        'CollectRhs'        , CollectRhs                                  , ...
        'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs                , ...
        'CoeffsParams'      , struct('a',1,'b',1,'sigma',[])              , ...
        'ScattererHandle'   , @Tools.Scatterer.PolarScatterer             , ...
        'ScattererParams'   , struct('r0',[], 'Stencil', Stencil) , ...
        'DiffOp'            , DiffOp                                      , ...
        'DiffOpParams'      , DiffOpParams                                , ...
        'SourceHandle'      , @Tools.Source.NavierStokesDescreteSrc       , ...
        ... %'SourceParams'     , SourceParams                                , ...
        'Extension'         , Extension                                   , ...
        'ExtensionParams'   , []                                          , ...
        'ExtensionPsi'      , ExtensionPsi                                );
    
    
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



function MakeExactMovie(Grid,tN, Np,Params,Ft,ht,Exact)
    
    Record = Tools.Common.SimpleRecorder(Grid,tN);
    
    Ex = zeros(size(Grid.R));
    Params.r=Grid.R(Np);
    for tn=0:tN
        t = tn*ht;
        Ex(Np) = Exact.Omega(Grid.Theta(Np),Params)*Ft(t);
        
        Record.Store(t,Ex);
    end
    
    save(['ExOmega' func2str(Ft) '.mat'],'Record');
    
end




function g=g1(ExParam,k,Fn,Grid,ScattererParams, Exact)
    g = zeros(Grid.Nx,Grid.Ny);
    
   
    
    S = Tools.Scatterer.PolarScatterer(Grid,ScattererParams);
    Msk = S.Mp;

    ExParam.r = Grid.R(Msk);
    theta =Grid.Theta(Msk); 
    
    g(Msk) = Exact.L_np1(theta,ExParam, k,Fn,1);
end

function gnp1 = TimeStepSource(theta,Params,gn,DescreteOmega,TimeFunc,ht,n,RN,UseConvTerm,Exact,D)
    %here we calc g^{n+1} from g^n
    
    % 'Src Derivatives line 135' On Analytical Scatterer
    % 'time step update line 108' Mp
        
    k = 2*RN/ht;
    
    Src = @(n) Exact.NSTimeSource(theta,Params,RN,TimeFunc,n*ht,UseConvTerm);
    
    gnp1 = 2*RN*D -(gn + 2*k*DescreteOmega +  RN*(Src(n+1) + Src(n)) );
end

