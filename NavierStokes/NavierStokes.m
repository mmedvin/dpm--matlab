function NavierStokes
    
    LinearSolverType = 0;%3;%
    CollectRhs = 1;
    KindOfConvergance = Tools.Enums.Convergance.Exact;

    CoeffsParams = struct('a',1,'b',1,'sigma',123);
    ExactChoice = NavierStokesExact.Exact2DiscreteSrc;
    [Exact,SourceHandle] = ExactChoice.Helper();
    
    BasisType = Tools.Enums.Basis.Fourier;
    
    test= -1;
    switch test
        case 1
            TstName = 'Dirichlet for Omega';           
        case 2
            TstName = 'Neumann for Omega';
        case 3
            TstName = 'Dirichlet for Omega, but exact Omega for Psi';
        case 4
            TstName = 'Neumann for Omega, but exact xi_psi';
        otherwise
            TstName = 'NS';
            
    end
    
    
    
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
        
        Extension = @Tools.Extensions.EBPolarLaplace5OrderExtension;
        
        Stencil=13;
        DiffOp = @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs;
        
    end
    
    ExParams.r0=R0;
    for   r_ = R0*0.9%[0.9,1,1.1]
        ExParams.r = r_;
        Eparam =@(r) struct('r0',ExParams.r0,'r',r);
                
     
        f   =@(phi) Exact.Omega(phi,ExParams);
        fn  =@(phi) Exact.DrDOmega(phi,ExParams);
        Basis = BasisType.Helper(f,fn);
        
        ExtensionParamsPsi.PsiBC = PsiBC(Eparam,Exact);
        
        ErrInfPre = 0; Err2Pre = 0;  biPre=0; b2Pre=0; pbiPre=0; pb2Pre=0; ErrpiPre=0; Errp2Pre=0;
        
        fprintf('Navier Stokes, %s\n',ExactChoice.toString());
        fprintf('NBss0=%d, NBss1=%d, LinearSolverType = %d , Order=%d ,k=%-2.2f,Scatterer at r=%-2.2f, %s, %s\n', ...
            Basis.NBss0, Basis.NBss1, LinearSolverType,Order,CoeffsParams.sigma,ExParams.r,KindOfConvergance.toString(),TstName);
        
        for n=1:5 %6 %run different grids
            tic
            %build grid
            p=3;%4;
            Nx=2.^(n+p)+1;
            Ny=2.^(n+p)+1;
            
            Grid = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
            SourceParams = ExParams;
            if ExactChoice.Subtype == Tools.Enums.Category.Discrete 
                SourceParams = ExParams;
                SourceParams.Src = Exact.NSSource(Grid.Theta,Eparam(Grid.R), CoeffsParams.sigma);
            end
            
            %k=1;%2/abs(Grid.x(1)-Grid.x(2));
            
            Setup  = struct( 'Basis'     , Basis, ...
                'Grid'              , Grid, ...
                'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                'CoeffsParams'      , CoeffsParams, ...
                'ScattererHandle'   , @Tools.Scatterer.PolarScatterer, ...
                'ScattererParams'   , struct('r0',ExParams.r, 'Stencil', Stencil), ...
                'CollectRhs'        , CollectRhs, ...
                'DiffOp'            , DiffOp, ...
                'DiffOpParams'      , struct(   'BC_y1',  0, 'BC_yn',  0,'BC_x1',0 , 'BC_xn',0, 'LinearSolverType', LinearSolverType, 'Order',Order), ...
                'SourceHandle'      , SourceHandle, ...
                'SourceParams'      , SourceParams, ...
                'Extension'         , Extension, ...
                'ExtensionParams'   , [], ...
                'ExtensionPsi'      , ExtensionPsi, ...
                'ExtensionParamsPsi', ExtensionParamsPsi ...
                );
            
            
            Prb =  Solvers.NavierStokesSolver(Setup);

            
            switch test
                case {1,3} % dirichlet problem
                    cn =[Basis.cn0; ( Prb.Q1 \ ( -Prb.Q0*Basis.cn0 - Prb.TrGF{1} - Prb.Qf{1}))];
                case {2,4} %Neumann
                    cn =[( Prb.Q0 \ ( -Prb.Q1*Basis.cn1 - Prb.TrGF{1} - Prb.Qf{1})); Basis.cn1];                                       
                otherwise % navier stockes
                    rhs = [-Prb.Qf{1}-Prb.TrGF{1} ; -Prb.Qpsi{1} - Prb.QPsif - Prb.TrGGF - Prb.TrGpsiGExf ];
                    cn = [ Prb.Q{1},Prb.Q{2} ; Prb.QpsiOmega{1} + Prb.GPO{1}, Prb.QpsiOmega{2} + Prb.GPO{2}]\rhs;
            end
            
            Oxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
            Oxi(Prb.GridGamma) = [Prb.W{1}(Prb.GridGamma,:),Prb.W{2}(Prb.GridGamma,:)]*cn  + Prb.Wf{1}(Prb.GridGamma);
            
            Pxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
            Pxi(Prb.GridGamma) = [Prb.WpsiOmega{1}(Prb.GridGamma,:),Prb.WpsiOmega{2}(Prb.GridGamma,:)]*cn + Prb.Wpsi{1}(Prb.GridGamma,:) + Prb.WPsif{1}(Prb.GridGamma,:);
            
            omega = Prb.P_Omega(Oxi);
                        
            if test==3 %exOmega
                ExOmega = zeros(size(Grid.R));
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.Scatterer.Np);
                ExOmega(Prb.Scatterer.Np) = Omega(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                p = Prb.P_OmegaPsi(ExOmega,Pxi);
                
            elseif test==4 %exPxi 
                ExParams2 = ExParams;
                ExParams2.r = Prb.Scatterer.r;
                
                exPxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                exPxi(Prb.GridGamma) = Psi(Prb.Scatterer.th,ExParams2);
                
                p = Prb.P_OmegaPsi(omega,exPxi);
                
            else
                p = Prb.P_OmegaPsi(omega,Pxi);%O,Or,Ott);
            end
            
            t1=toc;
            
            %------------------------------------------------------------------
            % Comparison
            %------------------------------------------------------------------
            
            
            
            if strcmpi(KindOfConvergance,'Grid')
                if n > 1
                    u1= spalloc(Nx,Ny,nnz(u0));
                    u1(Prb.Np) = u0(Prb.Np);
                    
                    tmp = omega(1:2:end,1:2:end)-u1(1:2:end,1:2:end);
                    
                    %etinf(n) =norm(tmp(:),inf);
                    ErrInf = norm(tmp(:),inf);
                    %fprintf('k=%d,M=%d,N=%-10dx%-10d\t etinf=%d\ttime=%d\n',k,Basis.M, Nr,Nth,full(etinf(n)),t);
                    %fprintf('k=%d,NBss0=%d,NBss1=%d,N=%-10dx%-10d\t ErrTot=%d\t rate=%-5.2f\t time=%d\n',k,Basis.NBss0,Basis.NBss1, Nr,Nth,ErrTot,log2(ErrPre/ErrTot),t);
                    fprintf('N=%-10dx%-10d\t ErrTot=%d\t rate=%-5.2f\t time=%d\n', Nx,Ny,ErrInf,log2(ErrInfPre/ErrInf),t1);
                    
                    ErrInfPre = ErrInf;
                    
                    
                    % fprintf('k=%-5.4f,pf=%-5.4f,pf=%-5.4f,IncAng=%d,M=%d,N=%-10dx%-10d\t etinf=%d\ttime=%d\t \n',k,  (k^5)* Grid.dx^4, (k^5) * Grid.dy^4,IncAngD,Basis.M, Nr,Nth,full(etinf(n)),t);
                    
                end
                
                u0=spalloc(Nx*2-1,Ny*2-1,nnz(omega));
                u0(1:2:end,1:2:end)=omega;
                
            elseif strcmpi(KindOfConvergance,'Exact')
                ExParams2 = ExParams;
                ExParams2.r = Prb.Scatterer.r;
                OxiEx = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                OxiEx(Prb.GridGamma) = Exact.Omega(Prb.Scatterer.th,ExParams2);
                [bi,b2] = cmpr(OxiEx(Prb.GridGamma), Oxi(Prb.GridGamma));
                
                PxiEx = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                PxiEx(Prb.GridGamma) = Exact.Psi(Prb.Scatterer.th,ExParams2);
                [pbi,pb2] = cmpr(PxiEx(Prb.GridGamma), Pxi(Prb.GridGamma));
                
                
                Ex = zeros(size(Grid.R));
                
                %Ex(IntPrb.Scatterer.Np) = Ex(FocalDistance,Prb.Scatterer.R(Prb.Scatterer.Np),Prb.Scatterer.Th(Prb.Scatterer.Np),R0);
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.Scatterer.Np);
                Ex(Prb.Scatterer.Np) = Exact.Omega(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                
                ExP = zeros(size(Grid.R));
                ExP(Prb.Scatterer.Np) = Exact.Psi(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                
                [ErrInf,Err2] = cmpr(Ex(Prb.Scatterer.Np),omega(Prb.Scatterer.Np));
                [Errpi,Errp2] = cmpr(ExP(Prb.Scatterer.Np),p(Prb.Scatterer.Np));
                
                str = sprintf('$%-6d\\\\times%-7d$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ \\\\\\\\ \n',...
                    Nx,Ny,ErrInf, log2(ErrInfPre/ErrInf), ...
                    Errpi,  log2(ErrpiPre/Errpi)  , ...
                    bi,  log2(biPre/bi) , ...
                    pbi,  log2(pbiPre/pbi)   )     ;
                
                ErrInfPre = ErrInf;
                Err2Pre = Err2;
                biPre=bi;
                b2Pre=b2;
                pbiPre=pbi;
                pb2Pre=pb2;
                ErrpiPre = Errpi;
                Errp2Pre = Errp2;
                
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
   % i=find(abs(tmp)==Linf,1);
end



