function NavierStokes2
    
    for test=-1%[3,4]%[-1,1,2];
        switch test
            case 1
                TstName = 'Dirichlet for Omega';
                
            case 2
                TstName = 'Neumann for Omega';
            case 3
                TstName = 'Dirichlet for Omega, but exact Omega for Psi';
            case 4
                TstName = 'Neumann for Omega, but exact Omega for Psi';
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
            DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
            %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
        elseif Order==4
            
            Extension = @Tools.Extensions.EBPolarLaplace5OrderExtension;
            
            Stencil=13;
            DiffOp = @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs;
            
        end
        
        ExParams.r0=R0;
        ExParams.r=R0*1;
        ExParams.p=99;
        Eparam =@(r) struct('r0',ExParams.r0,'r',r,'p',ExParams.p);
        
        
        k=10;
        
        BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
        ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
        
        f   =@(phi) Omega(phi,ExParams);
        fn  =@(phi) DrDOmega(phi,ExParams);
        
        xi0Psi  =@(theta,r) Psi(theta,Eparam(r));
        xi1Psi  =@(theta,r) DrDPsi(theta,Eparam(r));
        xi0PsiTT=@(theta,r) DPsiDThetaTheta(theta,Eparam(r));
        
        
        if strcmpi(BType,'Chebyshev')
            Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,fn,ChebyshevRange);
        elseif strcmpi(BType,'Fourier')
            Basis = Tools.Basis.FourierBasis.BasisHelper(f,fn);%,[3,3]);%,[15,15]);%[1e-14,1e-14]);%20);%
        end
        
        LinearSolverType = 0;
        CollectRhs = 1;
        
        ErrInfPre = 0; Err2Pre = 0;
        biPre=0; b2Pre=0;
        ErrpiPre=0; Errp2Pre=0;
        
        fprintf('Navier Stokes, NBss0=%d, NBss1=%d, LinearSolverType = %d , Order=%d ,k=%-2.2f,Scatterer at r=%-2.2f, %s, %s\n', ...
            Basis.NBss0, Basis.NBss1, LinearSolverType,Order,k,ExParams.r,strKoC,TstName);
        
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
                'SourceHandle'      , @Tools.Source.NavierStokesSourceRTh2, ...
                'SourceParams'      , ExParams, ...
                'Extension'         , Extension, ...
                'ExtensionParams'   , [], ...
                'PsiBC'             , struct('xi0Psi',xi0Psi,'xi1Psi',xi1Psi,'xi0PsiTT',xi0PsiTT)...
                );
            
            
            Prb =  Solvers.NavierStokesSolver(Setup);
            %InteriorLaplacianSolver(Setup);
            %;
            
            switch test
                case {1,3} % dirichlet problem
                    cn =[Basis.cn0; ( Prb.Q1 \ ( -Prb.Q0*Basis.cn0 - Prb.TrGF{1} - Prb.Qf{1}))];
                case {2,4} %Neumann
                    cn =[( Prb.Q0 \ ( -Prb.Q1*Basis.cn1 - Prb.TrGF{1} - Prb.Qf{1})); Basis.cn1];
                case 99
                    % cn =[Basis.cn0; ( Prb.Q1 \ ( -Prb.Q0*Basis.cn0 - Prb.TrGF{1} - Prb.Qf{1}))];
                    % obj.OpPsi.Solve(w-GLW);
                    % rhs = [-Prb.Qf{1}-Prb.TrGF{1} ;  -Prb.TrGGF ];
                    % cn = [ Prb.Q{1},Prb.Q{2} ; Prb.P{1},Prb.P{2}]\rhs;
                    
                otherwise % navier stockes
                    rhs = [-Prb.Qf{1}-Prb.TrGF{1} ; -Prb.TrGGF ];
                    cn = [ Prb.Q{1},Prb.Q{2} ; Prb.P{1},Prb.P{2}]\rhs;
            end
            if 1
                %rhs = [-Prb.Qf{1}-Prb.TrGF{1} + Prb.TrGGF; zeros(numel(Prb.GridGamma),1) ];
                
            else % debuging
                
            end
            
            Oxi = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
            Oxi(Prb.GridGamma) = [Prb.W{1}(Prb.GridGamma,:),Prb.W{2}(Prb.GridGamma,:)]*cn  + Prb.Wf{1}(Prb.GridGamma);
            
            
            u = Prb.P_Omega(Oxi);
            
            O=0;
            
            for j=1:numel(Basis.Indices0)
                uj = Basis.Handle(Prb.Scatterer.BasisArg,Basis.Indices0(j),[]);
                O = O + cn(j).*uj.xi0;
            end
            
            if test==3
                exu = zeros(size(Grid.R));
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.Scatterer.Np);
                exu(Prb.Scatterer.Np) = Omega(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                p = Prb.SPsi(exu,O);
            elseif test==4
                
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.GridGamma);
                exO = Omega(Prb.Scatterer.Th(Prb.GridGamma),ExParams3);
                p = Prb.SPsi(u,exO);
            else
                p = Prb.SPsi(u,O);
            end
            
            t1=toc;
            
            %------------------------------------------------------------------
            % Comparison
            %------------------------------------------------------------------
            
            
            
            if strcmpi(KindOfConvergance,'Grid')
                if n > 1
                    u1= spalloc(Nx,Ny,nnz(u0));
                    u1(Prb.Np) = u0(Prb.Np);
                    
                    tmp = u(1:2:end,1:2:end)-u1(1:2:end,1:2:end);
                    
                    %etinf(n) =norm(tmp(:),inf);
                    ErrInf = norm(tmp(:),inf);
                    %fprintf('k=%d,M=%d,N=%-10dx%-10d\t etinf=%d\ttime=%d\n',k,Basis.M, Nr,Nth,full(etinf(n)),t);
                    %fprintf('k=%d,NBss0=%d,NBss1=%d,N=%-10dx%-10d\t ErrTot=%d\t rate=%-5.2f\t time=%d\n',k,Basis.NBss0,Basis.NBss1, Nr,Nth,ErrTot,log2(ErrPre/ErrTot),t);
                    fprintf('N=%-10dx%-10d\t ErrTot=%d\t rate=%-5.2f\t time=%d\n', Nx,Ny,ErrInf,log2(ErrInfPre/ErrInf),t1);
                    
                    ErrInfPre = ErrInf;
                    
                    
                    % fprintf('k=%-5.4f,pf=%-5.4f,pf=%-5.4f,IncAng=%d,M=%d,N=%-10dx%-10d\t etinf=%d\ttime=%d\t \n',k,  (k^5)* Grid.dx^4, (k^5) * Grid.dy^4,IncAngD,Basis.M, Nr,Nth,full(etinf(n)),t);
                    
                end
                
                u0=spalloc(Nx*2-1,Ny*2-1,nnz(u));
                u0(1:2:end,1:2:end)=u;
                
            elseif strcmpi(KindOfConvergance,'Exact')
                ExParams2 = ExParams;
                ExParams2.r = Prb.Scatterer.r;
                OxiEx = spalloc(Nx,Ny   ,numel(Prb.GridGamma));
                OxiEx(Prb.GridGamma) = Omega(Prb.Scatterer.th,ExParams2);
                [bi,b2] = cmpr(OxiEx(Prb.GridGamma), Oxi(Prb.GridGamma));
                
                Ex = zeros(size(Grid.R));
                
                %Ex(IntPrb.Scatterer.Np) = Ex(FocalDistance,Prb.Scatterer.R(Prb.Scatterer.Np),Prb.Scatterer.Th(Prb.Scatterer.Np),R0);
                ExParams3 = ExParams;
                ExParams3.r = Prb.Scatterer.R(Prb.Scatterer.Np);
                Ex(Prb.Scatterer.Np) = Omega(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                
                ExP = zeros(size(Grid.R));
                ExP(Prb.Scatterer.Np) = Psi(Prb.Scatterer.Th(Prb.Scatterer.Np),ExParams3);
                
                [ErrInf,Err2] = cmpr(Ex(Prb.Scatterer.Np),u(Prb.Scatterer.Np));
                [Errpi,Errp2] = cmpr(ExP(Prb.Scatterer.Np),p(Prb.Scatterer.Np));
                
                str = sprintf('$%-6d\\\\times%-7d$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ \\\\\\\\ \n',...
                    Nx,Ny,ErrInf, log2(ErrInfPre/ErrInf), ...
                    Errpi,  log2(ErrpiPre/Errpi)  , ...
                    bi,  log2(biPre/bi)          );
                
                ErrInfPre = ErrInf;
                Err2Pre = Err2;
                biPre=bi;
                b2Pre=b2;
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
end

function P = Psi(theta,Params)
    r = Params.r;
    x = r .* cos(theta);
    y = r .* sin(theta);
    
    P = sin(x).*sin(y);
end

function DP = DrDPsi(theta,Params)
    r = Params.r;
    x = r .* cos(theta);
    y = r .* sin(theta);
    
    xr = cos(theta);
    yr = sin(theta);
    
    DP = cos(x).*sin(y).*xr + cos(y).*sin(x).*yr;
end

function DP = DPsiDThetaTheta(theta,Params)
    r = Params.r;
    x = r .* cos(theta);
    y = r .* sin(theta);
    
    DP = -(r.*(sin(x).*(r.*sin(y) + cos(y).*sin(theta)) + cos(x).*(cos(theta).*sin(y) + r.*cos(y).*sin(2*theta))));
end

function O = Omega(theta,Params)
    
    O = -2*Psi(theta,Params);
    
end

function DO = DrDOmega(theta,Params)
    DO = -2*DrDPsi(theta,Params);
end

