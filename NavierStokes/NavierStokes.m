function NavierStokes
    
    LinearSolverType = 0;%3;%
    CollectRhs = 1;
    
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
    
    KindOfConvergance='Exact';%'Grid';%
    
    if strcmpi(KindOfConvergance,'Grid')
        strKoC = 'grid convergence';
    elseif strcmpi(KindOfConvergance,'Exact')
        strKoC = 'convergance to exact solution';
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
        ExtensionPsi = @Tools.Extensions.NavierStokesPsi4rdOrderExtension;
        
        DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
        %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
    elseif Order==4
        
        Extension = @Tools.Extensions.EBPolarLaplace5OrderExtension;
        
        Stencil=13;
        DiffOp = @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs;
        
    end
    
    ExParams.r0=R0;
    for   r_ = R0*[0.9,1,1.1]
        ExParams.r = r_;
        ExParams.p=99;
        Eparam =@(r) struct('r0',ExParams.r0,'r',r,'p',ExParams.p);
        
        
        k=1000;
        
        BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
        ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
        
        f   =@(phi) Omega(phi,ExParams);
        fn  =@(phi) DrDOmega(phi,ExParams);
        
        xi0Psi  =@(theta,r) Psi(theta,Eparam(r));
        xi1Psi  =@(theta,r) DrDPsi(theta,Eparam(r));
        xi0PsiTT=@(theta,r) DPsiDThetaTheta(theta,Eparam(r));
        xi1PsiTT=@(theta,r) DPsiDrThetaTheta(theta,Eparam(r));
        xi0PsiTTTT=@(theta,r) DPsiD4Theta(theta,Eparam(r));
        
        ExtensionParamsPsi.PsiBC = struct('xi0Psi',xi0Psi,'xi1Psi',xi1Psi,'xi0PsiTT',xi0PsiTT,'xi1PsiTT',xi1PsiTT,'xi0PsiTTTT',xi0PsiTTTT);
        
        if strcmpi(BType,'Chebyshev')
            Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,fn,ChebyshevRange);
        elseif strcmpi(BType,'Fourier')
            Basis = Tools.Basis.FourierBasis.BasisHelper(f,fn);%,[3,3]);%,[15,15]);%[1e-14,1e-14]);%20);%
        end
        
        ErrInfPre = 0; Err2Pre = 0;  biPre=0; b2Pre=0; pbiPre=0; pb2Pre=0; ErrpiPre=0; Errp2Pre=0;
        
        fprintf('Navier Stokes, NBss0=%d, NBss1=%d, LinearSolverType = %d , Order=%d ,k=%-2.2f,Scatterer at r=%-2.2f, %s, %s\n', Basis.NBss0, Basis.NBss1, LinearSolverType,Order,k,ExParams.r,strKoC,TstName);
        
        for n=1:5 %6 %run different grids
            tic
            %build grid
            p=3;%4;
            Nx=2.^(n+p)+1;
            Ny=2.^(n+p)+1;
            
            Grid = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
            
            %k=1;%2/abs(Grid.x(1)-Grid.x(2));
            
            Setup  = struct( 'Basis'     , Basis, ...
                'Grid'              , Grid, ...
                'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                'CoeffsParams'      , struct('a',1,'b',1,'sigma',k), ...
                'ScattererHandle'   , @Tools.Scatterer.PolarScatterer, ...
                'ScattererParams'   , struct('r0',ExParams.r, 'Stencil', Stencil), ...
                'CollectRhs'        , CollectRhs, ...
                'DiffOp'            , DiffOp, ...
                'DiffOpParams'      , struct(   'BC_y1',  0, 'BC_yn',  0,'BC_x1',0 , 'BC_xn',0, 'LinearSolverType', LinearSolverType, 'Order',Order), ...
                'SourceHandle'      , @Tools.Source.NavierStokesSourceRTh, ...
                'SourceParams'      , ExParams, ...
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
                    rhs = [-Prb.Qf{1}-Prb.TrGF{1} ; -Prb.Qpsi{1}-Prb.TrGGF - Prb.TrGpsiGExf ];
                    cn = [ Prb.Q{1},Prb.Q{2} ; Prb.QpsiOmega{1} + Prb.GPO{1}, Prb.QpsiOmega{2} + Prb.GPO{2}]\rhs;
            end
            
            Oxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
            Oxi(Prb.GridGamma) = [Prb.W{1}(Prb.GridGamma,:),Prb.W{2}(Prb.GridGamma,:)]*cn  + Prb.Wf{1}(Prb.GridGamma);
            
            Pxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
            Pxi(Prb.GridGamma) = [Prb.WpsiOmega{1}(Prb.GridGamma,:),Prb.WpsiOmega{2}(Prb.GridGamma,:)]*cn + Prb.Wpsi{1}(Prb.GridGamma,:);
            
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
                OxiEx(Prb.GridGamma) = Omega(Prb.Scatterer.th,ExParams2);
                [bi,b2] = cmpr(OxiEx(Prb.GridGamma), Oxi(Prb.GridGamma));
                
                PxiEx = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                PxiEx(Prb.GridGamma) = Psi(Prb.Scatterer.th,ExParams2);
                [pbi,pb2] = cmpr(PxiEx(Prb.GridGamma), Pxi(Prb.GridGamma));
                
                
                Ex = zeros(size(Grid.R));
                
                %Ex(IntPrb.Scatterer.Np) = Ex(FocalDistance,Prb.Scatterer.R(Prb.Scatterer.Np),Prb.Scatterer.Th(Prb.Scatterer.Np),R0);
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.Scatterer.Np);
                Ex(Prb.Scatterer.Np) = Omega(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                
                ExP = zeros(size(Grid.R));
                ExP(Prb.Scatterer.Np) = Psi(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                
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
