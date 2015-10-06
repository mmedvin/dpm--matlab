                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          function RunLaplacian351
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    a=1;%2.5;
    b=1/2;	

	
	%y1=-0.7;yn=0.7;
	%y1=-1.1;yn=1.1;
    
    %x1i=-1.2;xni=1.2;
    %y1i=-0.7;yni=0.7; Dom='Diff';

    %x1i=-1.1;xni=1.1;
    %y1i=-0.6;yni=0.6; Dom='Diff';

    
    x1i=-2;xni=2;
     y1i=-2;yni=2; Dom='Same';
  
    %x1e=-1.2;xne=1.2;
	%y1e=-1.2;yne=1.2;
    
    x1e=-2;xne=2;
	y1e=-2;yne=2;

    BIn = 1e0;
    BOut = 1e3;

	BType		= 'Fourier';
    
    for Order= [2,4]
        if Order==2 , ExpansionType=33; Stencil=5; else ExpansionType=35; Stencil=9; end
                
        for b=1/2%[1/2,1/16] %[1/2,1/4,1/8,1/16]
            
            %y1i=-(b+0.2);yni=(b+0.2); Dom='Diff';
            
            FocalDistance = sqrt(a^2-b^2);
            Eta0 = acosh(a/FocalDistance);
            
            f   =@(phi) ExtExact(FocalDistance,Eta0,phi);
            g   =@(phi) IntExact(FocalDistance,Eta0,phi);
            gn   =@(phi) dndIntExact(FocalDistance,Eta0,phi);
            fn   =@(phi) dndExtExact(FocalDistance,Eta0,phi);
                                   
            if strcmpi(BType,'Chebyshev')
                Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
            elseif strcmpi(BType,'Fourier')
                Basis = Tools.Basis.FourierBasis.BasisHelper(@(phi) g(phi)-f(phi),@(phi) BIn*gn(phi) - BOut*fn(phi));
            end
            
            for RelErr=[0,1]
                
                for	   LinearSolverType = 0
                    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end
                    
                    % %         ErrIntPre = 0; 	ErrExtPre = 0;	%ErrTotPre = 0; ErrUPre = 0; ErrUrPre = 0; ErrUrrPre = 0;
                    % %         ErrUPre = 0; ErrUxPre = 0; ErrUxxPre = 0;  ErrUyPre = 0; ErrUyyPre = 0;  ErrUxyPre = 0;
                    
                    ErrUInfPre = 0; ErrU2Pre = 0; ErrUxInfPre = 0; ErrUx2Pre = 0; ErrUyInfPre = 0; ErrUy2Pre = 0; ErrUxxInfPre = 0; ErrUxx2Pre = 0; ErrUyyInfPre = 0; ErrUyy2Pre = 0; ErrUxyInfPre = 0; ErrUxy2Pre = 0;
                    
                    PreAmble = sprintf('%%%%Problem 3.51, RelErr=%d, NBss0=%d,, NBss1=%d Order = %d, AR=%d, LinearSolverType = %d, BIn=%d, BOut=%d, ID=[%d,%d]x[%d,%d], ED=[%d,%d]x[%d,%d]\n',...
                        RelErr, Basis.NBss0, Basis.NBss1, Order,a/b, LinearSolverType, BIn, BOut,x1i,xni,y1i,yni,x1e,xne,y1e,yne);
                    fprintf(PreAmble);
                    
                    %FileName = sprintf('P351RE%dNBss0%dNBss1%dO%dAR%dSol%dBI%dBO%d%sDom',RelErr, Basis.NBss0, Basis.NBss1, Order,a/b, LinearSolverType, BIn, BOut,Dom);
                    %fileID = fopen([FileName '.txt'],'a');
                    %fprintf(fileID,PreAmble);
                    
                    
                    for n=1:6 %run different grids
                        tic
                        %build grid
                        % 		p=3;%4;
                        % 		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
                        
                        p=20;
                        Nx=p*2^n; Ny=p*2^n;
                        
                        GridInt                = Tools.Grid.CartesianGrid(x1i,xni,Nx,y1i,yni,Ny);
                        
                        ScattererHandle  = @Tools.Scatterer.EllipticScatterer;
                        ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',ExpansionType, 'Stencil', Stencil);
                        
                        %------------------------------------------------------------------
                        
                        if Order==4
                            DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat4OrdrConst;%LaplacianOpBCinMat4OrdrConst; %LaplacianOpBCinMat;
                        elseif Order==2
                            DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
                            %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
                        end
                        
                        IntPrb = Solvers.InteriorHomoLaplacianSolver(struct(...
                            'Basis'             , Basis, ...
                            'Grid'              , GridInt, ...
                            'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                            'CoeffsParams'      , struct('a',BIn,'b',BIn,'sigma',0), ...
                            'ScattererHandle'   , ScattererHandle, ...
                            'ScattererParams'   , ScattererParams, ...
                            'CollectRhs'        , CollectRhs, ...
                            'DiffOp'            , DiffOp, ...
                            'DiffOpParams'      , struct('BC_x1',0,'BC_xn', 0,'BC_y1',0,'BC_yn',0, 'LinearSolverType', LinearSolverType), ...
                            'Extension'         , @Tools.Extensions.TwoTupleExtension, ...
                            'ExtensionParams',[] ...
                            ));
                        
                        %------------------------------------------------------------------
                        GridExt                = Tools.Grid.CartesianGrid(x1e,xne,Nx,y1e,yne,Ny);
                        
                        if isequal(DiffOp , @Tools.DifferentialOps.LaplacianOpBCinRhs)
                            DiffOpParamsExt = struct('BC_y1', sin(GridExt.x(1)).*cos(GridExt.y(2:end-1)).','BC_yn',  sin(GridExt.x(end)).*cos(GridExt.y(2:end-1)).',...
                                'BC_x1',(sin(GridExt.x(2:end-1)).*cos(GridExt.y(1))),'BC_xn',(sin(GridExt.x(2:end-1)).*cos(GridExt.y(end))), 'LinearSolverType', LinearSolverType);
                        elseif isequal(DiffOp , @Tools.DifferentialOps.LaplacianOpBCinMat) || isequal(DiffOp , @Tools.DifferentialOps.LaplacianOpBCinMat4OrdrConst)
                            DiffOpParamsExt = struct('BC_y1', sin(GridExt.x(1)).*cos(GridExt.y).','BC_yn',  sin(GridExt.x(end)).*cos(GridExt.y).',...
                                'BC_x1',(sin(GridExt.x).*cos(GridExt.y(1))),'BC_xn',(sin(GridExt.x).*cos(GridExt.y(end))), 'LinearSolverType', LinearSolverType);
                        end
                        
                        %ExtPrb =  Solvers.ExtLapSolver4Ordr( struct(...
                        ExtPrb =  Solvers.ExteriorLaplacianSolver( struct(...
                            'Basis'             , Basis, ...
                            'Grid'              , GridExt, ...
                            'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                            'CoeffsParams'      , struct('a',BOut,'b',BOut,'sigma',0), ...
                            'ScattererHandle'   , ScattererHandle, ...
                            'ScattererParams'   , ScattererParams, ...
                            'CollectRhs'        , CollectRhs, ...
                            'DiffOp'            , DiffOp, ...
                            'DiffOpParams'      , DiffOpParamsExt, ...
                            'SourceHandle'      , @Tools.Source.LaplaceSource_IIM351_Exterior, ...
                            'SourceParams'      , [], ...
                            'Extension'         , @Tools.Extensions.TwoTupleExtension, ...
                            'ExtensionParams',[] ...
                            ));
                        
                        %------------------------------------------------------------------
                        if 1
                            nGGext = numel(ExtPrb.GridGamma);
                            nGGint = numel(IntPrb.GridGamma);
                            
                            ZerosInt0=spalloc(nGGint,Basis.NBss0,0);
                            ZerosInt1=spalloc(nGGint,Basis.NBss1,0);
                            ZerosExt0=spalloc(nGGext,Basis.NBss0,0);
                            ZerosExt1=spalloc(nGGext,Basis.NBss1,0);
                            
                            Eye0 = speye(Basis.NBss0);
                            Eye1 = speye(Basis.NBss1);
                            Zeros2_0=spalloc(Basis.NBss1,Basis.NBss0,0);
                            Zeros2_1=spalloc(Basis.NBss0,Basis.NBss1,0);
                            
                            
                            rhs = zeros(nGGext+nGGint + Basis.NBss0+Basis.NBss1,1);
                            %rhs(1:nGG)	= (-IntPrb.TrGF -IntPrb.Qf);
                            rhs((nGGint+1):(nGGext+nGGint))= -(ExtPrb.TrGF +ExtPrb.Qf);
                            
                            x0 = @(phi) FocalDistance*cosh(Eta0).*cos(phi);
                            y0 = @(phi) FocalDistance*sinh(Eta0).*sin(phi);
                            
                            xn0 = @(phi) FocalDistance*sinh(Eta0).*cos(phi);
                            yn0 = @(phi) FocalDistance*cosh(Eta0).*sin(phi);
                            
                            ICu = @(phi) g(phi) - f(phi);%(x0(phi).^2 - y0(phi).^2) - sin(x0(phi)).* cos(y0(phi));%@(phi) f(phi)-g(phi); %
                            ICuc = Tools.Basis.FourierBasis.Coeffs(ICu, Basis.Core,0);
                            
                            ICun = @(phi)  BIn*gn(phi) - BOut*fn(phi);%2*BIn*(x0(phi).*xn0(phi)-y0(phi).*yn0(phi)) - BOut*(cos(x0(phi)).* cos(y0(phi)).*xn0(phi) - sin(x0(phi)).* sin(y0(phi)).*yn0(phi)); @(phi) BIn*fn(phi)-BOut*gn(phi);%
                            %2*BOut* sin(x0).* cos(y0);
                            ICunc = Tools.Basis.FourierBasis.Coeffs(ICun, Basis.Core,1);
                            
                            rhs((nGGext+nGGint+1):(nGGext+nGGint+Basis.NBss0))= ICuc;
                            rhs((nGGext+nGGint+Basis.NBss0+1):end)= ICunc;
                            
                            cn = [IntPrb.Q0,IntPrb.Q1,ZerosInt0, ZerosInt1;
                                ZerosExt0, ZerosExt1,ExtPrb.Q0,ExtPrb.Q1;
                                Eye0     , Zeros2_1 ,-Eye0     , Zeros2_1;
                                ...Zeros2_0,  Eye1, Zeros2_0, -Eye1]\rhs;
                                Zeros2_0, BIn*Eye1, Zeros2_0, -BOut*Eye1]\rhs;
                            
                            Intcn=cn(1:(Basis.NBss0+Basis.NBss1));
                            Extcn= cn((Basis.NBss0+Basis.NBss1+1):end);
                            
                            Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
                            Intxi(IntPrb.GridGamma) = [IntPrb.W{1}.W(IntPrb.GridGamma,:),IntPrb.W{2}.W(IntPrb.GridGamma,:)]*Intcn;% + IntPrb.Wf.W(IntPrb.GridGamma);
                            Intu = IntPrb.P_Omega(Intxi);
                            
                            Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
                            Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}.W(ExtPrb.GridGamma,:),ExtPrb.W{2}.W(ExtPrb.GridGamma,:)]*Extcn + ExtPrb.Wf.W(ExtPrb.GridGamma);
                            Extu =  spalloc(Nx,Ny,length(ExtPrb.GridGamma));
                            tmp = ExtPrb.P_Omega(Extxi);
                            Extu(ExtPrb.Scatterer.Nm) = tmp(ExtPrb.Scatterer.Nm);
                        else
                            %exterior
                            
                            EQ0 = ExtPrb.Q0;
                            EQ1 = ExtPrb.Q1;
                            ETrGF = ExtPrb.TrGF;
                            EQf = ExtPrb.Qf;
                            Extcn1 =( EQ1 \ ( -EQ0*Basis.cn0 -ETrGF - EQf)) ;
                            
                            Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
                            Extxi(ExtPrb.GridGamma) = ...
                                ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*Extcn1 + ExtPrb.Wf.W(ExtPrb.GridGamma);
                            
                            xiex  = ExtExact(FocalDistance,ExtPrb.Scatterer.eta,ExtPrb.Scatterer.phi);
                            Extxi(ExtPrb.GridGamma) = xiex; %debig
                            
                            Extu = ExtPrb.P_Omega(Extxi);
                            
                            %interior
                            
                            Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
                            %xiex = spalloc(Nx,Ny,length(IntPrb.GridGamma));
                            
                            % Neumann problem - not working need to add data for uniquiness
                            %Intcn0 =( IntPrb.Q0 \ ( -IntPrb.Q1*Basis.cn1 )) ;
                            %Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Intcn0 + IntPrb.W1(IntPrb.GridGamma,:)*Basis.cn1;
                            %xiex(IntPrb.GridGamma) = IntExact(FocalDistance,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi);
                            
                            Intcn0 = Basis.cn1;
                            IQ0 = IntPrb.Q0;
                            IQ1 = IntPrb.Q1;
                            Intcn1 =( IQ1 \ ( -IQ0*Intcn0 )) ;
                            Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Intcn0 + IntPrb.W1(IntPrb.GridGamma,:)*Intcn1;
                            
                            %xiex  = IntExact(FocalDistance,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi);
                            %Intxi(ExtPrb.GridGamma) = xiex; %debig
                            
                            Intu = IntPrb.P_Omega(Intxi);
                        end
                        % 		u=zeros(size(Grid.R));
                        % 		u(IntPrb.Np) = Intu(IntPrb.Np);
                        % 		u(ExtPrb.Nm) = Extu(ExtPrb.Nm);
                        
                        % 		u(IntPrb.Scatterer.Mm) = Intu(IntPrb.Scatterer.Mm);
                        % 		u(ExtPrb.Scatterer.Mp) = Extu(ExtPrb.Scatterer.Mp);
                        
                        
                        %	u(IntPrb.GridGamma) = (Intu(IntPrb.GridGamma) + Extu(IntPrb.GridGamma))/2;
                        %------------------------------------------------------------------
                        
                        t1=toc;
                        
                        %------------------------------------------------------------------
                        % Comparison
                        %------------------------------------------------------------------
                        
                        Intexact = zeros(size(GridInt.R));
                        Extexact = zeros(size(GridExt.R));
                        
                        Intexact(IntPrb.Scatterer.Np) = IntExact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np));
                        Extexact(ExtPrb.Scatterer.Nm) = ExtExact(FocalDistance,ExtPrb.Scatterer.Eta(ExtPrb.Scatterer.Nm),ExtPrb.Scatterer.Phi(ExtPrb.Scatterer.Nm));
                        
                        
                        Extu(1,:)   = Extexact(1,:);
                        Extu(end,:) = Extexact(end,:);
                        Extu(:,1)   = Extexact(:,1);
                        Extu(:,end) = Extexact(:,end);
                        
                        Iu=Intu;
                        %Iu(IntPrb.GridGamma)=0;
                        Ie=Intexact;
                        %Ie(IntPrb.GridGamma)=0;
                        
                        Eu=Extu;
                        %Eu(ExtPrb.GridGamma)=0;
                        Ee=Extexact;
                        %Ee(ExtPrb.GridGamma)=0;
                        
                        [ErrUInf,ErrU2] = cmpr([Iu(IntPrb.Scatterer.Mp);Eu(ExtPrb.Scatterer.Mm)],[Ie(IntPrb.Scatterer.Mp);Ee(ExtPrb.Scatterer.Mm)],[],RelErr);
                        
                        CDExt = Tools.Common.SecondDerivative(GridExt.Nx,GridExt.Ny,GridExt.dx,GridExt.dy);
                        CDInt = Tools.Common.SecondDerivative(GridInt.Nx,GridInt.Ny,GridInt.dx,GridInt.dy);
                        
                        [eTux,eTuy,eTuxx,eTuyy,eTuxy]       = CDExt.CartesianDerivatives(Extu);
                        [iTux,iTuy,iTuxx,iTuyy,iTuxy]       = CDInt.CartesianDerivatives(Intu);
                        [eExux,eExuy,eExuxx,eExuyy,eExuxy]  = CDExt.CartesianDerivatives(Extexact);
                        [iExux,iExuy,iExuxx,iExuyy,iExuxy]  = CDInt.CartesianDerivatives(Intexact);
                        
                        [ErrUxInf,ErrUx2]   = cmpr([iTux(IntPrb.Scatterer.Mp);eTux(ExtPrb.Scatterer.Mm)],[iExux(IntPrb.Scatterer.Mp);eExux(ExtPrb.Scatterer.Mm)],[],RelErr);
                        [ErrUyInf,ErrUy2]   = cmpr([iTuy(IntPrb.Scatterer.Mp);eTuy(ExtPrb.Scatterer.Mm)],[iExuy(IntPrb.Scatterer.Mp);eExuy(ExtPrb.Scatterer.Mm)],[],RelErr);
                        
                        [ErrUxxInf,ErrUxx2] = cmpr([iTuxx(IntPrb.Scatterer.Mp);eTuxx(ExtPrb.Scatterer.Mm)],[iExuxx(IntPrb.Scatterer.Mp);eExuxx(ExtPrb.Scatterer.Mm)],[],RelErr);
                        [ErrUyyInf,ErrUyy2] = cmpr([iTuyy(IntPrb.Scatterer.Mp);eTuyy(ExtPrb.Scatterer.Mm)],[iExuyy(IntPrb.Scatterer.Mp);eExuyy(ExtPrb.Scatterer.Mm)],[],RelErr);
                        [ErrUxyInf,ErrUxy2] = cmpr([iTuxy(IntPrb.Scatterer.Mp);eTuxy(ExtPrb.Scatterer.Mm)],[iExuxy(IntPrb.Scatterer.Mp);eExuxy(ExtPrb.Scatterer.Mm)],[],RelErr);
                        %------------------------------------------------------------------
                        %
                        %     fprintf('N=%-6dx%-7d Eu=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Eux=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euxx=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euyy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euxy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f timeA=%-6.2f\n',...
                        %         Nx,Ny,ErrUInf,   ErrU2,   log2(ErrUInfPre/ErrUInf),     log2(ErrU2Pre/ErrU2), ...
                        %               ErrUxInf,  ErrUx2,  log2(ErrUxInfPre/ErrUxInf),   log2(ErrUx2Pre/ErrUx2), ...
                        %               ErrUyInf,  ErrUy2,  log2(ErrUyInfPre/ErrUyInf),   log2(ErrUy2Pre/ErrUy2), ...
                        %               ErrUxxInf, ErrUxx2, log2(ErrUxxInfPre/ErrUxxInf), log2(ErrUxx2Pre/ErrUxx2), ...
                        %               ErrUyyInf, ErrUyy2, log2(ErrUyyInfPre/ErrUyyInf), log2(ErrUyy2Pre/ErrUyy2), ...
                        %               ErrUxyInf, ErrUxy2, log2(ErrUxyInfPre/ErrUxyInf), log2(ErrUxy2Pre/ErrUxy2), ...
                        %             t1);
                        
                        %$16\times16$	&	$9.8821e-04$&	$-	 $ 	&	$3.3575e-04$&	$-	 $	&	$2.7225e-03$&	$-	 $ 	&	$3.6927e-05$&	$-	 $ 	&	$2.3633e-05$&	$-	 $	&	$4.6294e-05$&	$-	 $	\\
                        
                        %  str = sprintf('$%-6d\\\\times%-7d$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ \\\\\\\\ \n',...
                        %        Nx-1,Ny-1,ErrUInf,  log2(ErrUInfPre/ErrUInf), ErrUxInf,log2(ErrUxInfPre/ErrUxInf), ErrUyInf,log2(ErrUyInfPre/ErrUyInf),...
                        %                    ErrU2,  log2(ErrU2Pre/ErrU2),     ErrUx2,   log2(ErrUx2Pre/ErrUx2), ErrUy2,  log2(ErrUy2Pre/ErrUy2) );
                        
                        str = sprintf('$%-6d\\\\times%-7d$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$  \\\\\\\\ \n',...
                            Nx,Ny,ErrUInf,  log2(ErrUInfPre/ErrUInf), ErrUxInf,log2(ErrUxInfPre/ErrUxInf), ErrUyInf,log2(ErrUyInfPre/ErrUyInf) );
                        % Nx-1,Ny-1,ErrUInf,  log2(ErrUInfPre/ErrUInf), ErrUxInf,log2(ErrUxInfPre/ErrUxInf), ErrUyInf,log2(ErrUyInfPre/ErrUyInf) );
                        
                        
                        fprintf(str);
                        %fprintf(fileID,str);
                        
                        
                        ErrUInfPre = ErrUInf; ErrU2Pre = ErrU2; ErrUxInfPre = ErrUxInf; ErrUx2Pre = ErrUx2; ErrUyInfPre = ErrUyInf; ErrUy2Pre = ErrUy2;
                        ErrUxxInfPre = ErrUxxInf; ErrUxx2Pre = ErrUxx2; ErrUyyInfPre = ErrUyyInf; ErrUyy2Pre = ErrUyy2;
                        ErrUxyInfPre = ErrUxyInf; ErrUxy2Pre = ErrUxy2;
                        
                        
                        %--------------------------
                        % %         	Ex=Tools.Exact.ExLapElps351(IntPrb.Scatterer, InteriorCoeffsHandle);
                        % %
                        % % 	u=Ex.u; %zeros(size(Grid.R));
                        % % 	%u(ExtPrb.Scatterer.Mm) = Extu(ExtPrb.Scatterer.Mm);
                        % % 	u(2:end-1,2:end-1) = Extu(2:end-1,2:end-1);
                        % % 	u(IntPrb.Scatterer.Mp) = Intu(IntPrb.Scatterer.Mp);
                        % % 	ErrU =  cmpr(Ex.u,u); %norm(Ex.u(:) - u(:),inf);
                        %--------------------------
                        % %     Intexact = zeros(size(GridInt.R));
                        % %     Extexact = zeros(size(GridExt.R));
                        % %
                        % %     Intexact(IntPrb.Scatterer.Np) = IntExact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np));
                        % %     Extexact(ExtPrb.Scatterer.Nm) = ExtExact(FocalDistance,ExtPrb.Scatterer.Eta(ExtPrb.Scatterer.Nm),ExtPrb.Scatterer.Phi(ExtPrb.Scatterer.Nm));
                        % %
                        % %      ErrInt = cmpr(Intexact(IntPrb.Scatterer.Np),Intu(IntPrb.Scatterer.Np));
                        % %      ErrExt = cmpr(Extexact(2:end-1,2:end-1),Extu(2:end-1,2:end-1));
                        
                        %--------------------------
                        
                        % %     CD = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy);
                        % % 	[Tux,Tuy,Tuxx,Tuyy,Tuxy] = CD.CartesianDerivatives(u(:));
                        % %
                        % %     Tux  = reshape(Tux,Grid.Nx,Grid.Ny);
                        % %     Tuy  = reshape(Tuy,Grid.Nx,Grid.Ny);
                        % %     Tuxx = reshape(Tuxx,Grid.Nx,Grid.Ny);
                        % %     Tuyy = reshape(Tuyy,Grid.Nx,Grid.Ny);
                        % %     Tuxy = reshape(Tuxy,Grid.Nx,Grid.Ny);
                        % %
                        % %     [Exux,Exuy,Exuxx,Exuyy,Exuxy] = CD.CartesianDerivatives(Ex.u(:));
                        % %
                        % %     Exux  = reshape(Exux,Grid.Nx,Grid.Ny);
                        % %     Exuy  = reshape(Exuy,Grid.Nx,Grid.Ny);
                        % %     Exuxx = reshape(Exuxx,Grid.Nx,Grid.Ny);
                        % %     Exuyy = reshape(Exuyy,Grid.Nx,Grid.Ny);
                        % %     Exuxy = reshape(Exuxy,Grid.Nx,Grid.Ny);
                        % %
                        % %     ux = Exux;
                        % %     ux(2:end-1,2:end-1) = Tux(2:end-1,2:end-1);
                        % %     ErrUx = cmpr(Exux,ux);%,ExtPrb.Scatterer.GridGamma);
                        % %
                        % %     uy = Exuy;
                        % %     uy(2:end-1,2:end-1) = Tuy(2:end-1,2:end-1);
                        % %     ErrUy = cmpr(Exuy,uy);%,ExtPrb.Scatterer.GridGamma);
                        % %
                        % %     uxx = Exuxx;
                        % %     uxx(2:end-1,2:end-1) = Tuxx(2:end-1,2:end-1);
                        % %     ErrUxx = cmpr(Exuxx,uxx);%,ExtPrb.Scatterer.GridGamma);
                        % %
                        % %     uyy = Exuyy;
                        % %     uyy(2:end-1,2:end-1) = Tuyy(2:end-1,2:end-1);
                        % %     ErrUyy = cmpr(Exuyy,uyy);%,ExtPrb.Scatterer.GridGamma);
                        % %
                        % %     uxy = Exuxy;
                        % %     uxy(2:end-1,2:end-1) = Tuxy(2:end-1,2:end-1);
                        % %     ErrUxy = cmpr(Exuxy,uxy);%,ExtPrb.Scatterer.GridGamma);
                        % %
                        % %     fprintf('N=%-6dx%-7d ErrInt=%-10.4d rt=%-6.2f ErrExt=%-10.4d rt=%-6.2f  Eu=%-10.4d rt=%-6.2f EUx=%d rt=%-6.2f EUy=%d rt=%-6.2f EUxx=%d rt=%-6.2f EUyy=%d rt=%-6.2f EUxy=%d rt=%-6.2f timeA=%-6.2f\n',...
                        % %         Nx,Ny,ErrInt,log2(ErrIntPre/ErrInt),ErrExt,log2(ErrExtPre/ErrExt),ErrU,log2(ErrUPre/ErrU),ErrUx,log2(ErrUxPre/ErrUx),ErrUy,log2(ErrUyPre/ErrUy),ErrUxx,log2(ErrUxxPre/ErrUxx),ErrUyy,log2(ErrUyyPre/ErrUyy),ErrUxy,log2(ErrUxyPre/ErrUxy),t1);
                        % %
                        % %     ErrUPre = ErrU; ErrUxPre = ErrUx; ErrUxxPre = ErrUxx;  ErrUyPre = ErrUy; ErrUyyPre = ErrUyy; ErrUxyPre = ErrUxy;
                        
                        %      dxdn = IntPrb.Scatterer.FocalDistance*sinh(IntPrb.Scatterer.Eta).*cos(IntPrb.Scatterer.Phi);
                        %      dydn = IntPrb.Scatterer.FocalDistance*cosh(IntPrb.Scatterer.Eta).*sin(IntPrb.Scatterer.Phi);
                        %
                        %     RD = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy,dxdn(:),dydn(:));
                        % 	[Tdudeta,Td2udeta2] = RD.RadialDerivatives(u(:));
                        %
                        %     Tdudeta  = reshape(Tdudeta,Grid.Nx,Grid.Ny);
                        %     Td2udeta2 = reshape(Td2udeta2,Grid.Nx,Grid.Ny);
                        %
                        %     [Exdudeta,Exd2udeta2] = RD.RadialDerivatives(Ex.u(:));
                        %
                        %     Exdudeta  = reshape(Exdudeta,Grid.Nx,Grid.Ny);
                        %     Exd2udeta2 = reshape(Exd2udeta2,Grid.Nx,Grid.Ny);
                        %
                        %
                        % 	dudeta = Exdudeta;%Ex.dudeta;
                        %     dudeta(2:end-1,2:end-1) = Tdudeta(2:end-1,2:end-1);
                        %     tmp = Exdudeta - dudeta;
                        % 	tmp(IntPrb.Scatterer.GridGamma) = 0;
                        % 	ErrUr = norm(tmp(:),inf);
                        %
                        % 	d2udeta2 = Exd2udeta2;%Ex.d2udeta2;
                        %     d2udeta2(2:end-1,2:end-1) = Td2udeta2(2:end-1,2:end-1);
                        %
                        % 	%tmp = Ex.d2udeta2 - d2udeta2;
                        %     tmp = Exd2udeta2 - d2udeta2;
                        % 	tmp(IntPrb.Scatterer.GridGamma) = 0;
                        %
                        % 	%tmp( IntPrb.Scatterer.Eta > (Eta0-0.2)  &  IntPrb.Scatterer.Eta < (Eta0+0.2) )=0;
                        % 	ErrUrr  = norm(tmp( :),inf);
                        % %------------------------------------------------------------------
                        %
                        % 	 Intexact = zeros(size(Grid.R));
                        % 	 Extexact = zeros(size(Grid.R));
                        %
                        %
                        % 	Intexact(IntPrb.Scatterer.Np) = IntExact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np));
                        % 	Extexact(ExtPrb.Scatterer.Nm) = ExtExact(FocalDistance,ExtPrb.Scatterer.Eta(ExtPrb.Scatterer.Nm),ExtPrb.Scatterer.Phi(ExtPrb.Scatterer.Nm));
                        %
                        %     ErrInt =norm(Intexact(IntPrb.Scatterer.Np)-Intu(IntPrb.Scatterer.Np),inf);
                        % 	%Extetinf(n) =norm(Extexact(ExtPrb.Scatterer.Nm)-Extu(ExtPrb.Scatterer.Nm),inf);
                        % 	tmp = Extexact(2:end-1,2:end-1)-Extu(2:end-1,2:end-1);
                        % 	ErrExt =norm(tmp(:),inf);
                        %
                        %     ErrTot = max(ErrInt,ErrExt);
                        %
                        %     fprintf('N=%-6dx%-7d EInt=%-10.4d rt=%-6.2f EExt=%-10.4d rt=%-6.2f ETot=%d\t rt=%-6.2f EU=%d\t rt=%-6.2f EUr=%d rt=%-6.2f EUrr=%d rt=%-6.2f  timeA=%-6.2f\n',...
                        %         Nx,Ny,ErrInt,log2(ErrIntPre/ErrInt),ErrExt,log2(ErrExtPre/ErrExt),ErrTot,log2(ErrTotPre/ErrTot),ErrU,log2(ErrUPre/ErrU),ErrUr,log2(ErrUrPre/ErrUr),ErrUrr,log2(ErrUrrPre/ErrUrr),t1);
                        %
                        % %          ErrIntPre = ErrInt;
                        % %          ErrExtPre = ErrExt;
                        %     ErrTotPre = ErrTot;
                        %     ErrUPre	  = ErrU;
                        %     ErrUrPre  = ErrUr;
                        %     ErrUrrPre  = ErrUrr;
                        
                    end
                end
            end
        end
    end
end



function [Linf,L2] = cmpr(ex,u,GG,RelErr)

    if nargin==2, GG=[];end
    tmp = ex - u;
    
    if ~isempty(GG)
        tmp(GG) = 0;
        u(GG)=0;
    end

    Linf = norm(tmp(:),inf);%/norm(u(:),inf);
    L2   = norm(tmp(:),2);%/norm(u(:),2);
    
    if RelErr
        Linf = Linf/norm(u(:),inf);
        L2   = L2/norm(u(:),2);
    end
end
	
function e = IntExact(FocalDist,eta,phi)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);
	
	e = x.^2 - y.^2;	
end
function dnde = dndIntExact(FocalDist,eta,phi)
	x  = FocalDist*cosh(eta).*cos(phi);
	y  = FocalDist*sinh(eta).*sin(phi);
	xn = FocalDist*sinh(eta).*cos(phi);
	yn = FocalDist*cosh(eta).*sin(phi);
	
	dnde = 2*x.*xn - 2.*y.*yn;
end

function e = ExtExact(FocalDist,eta,phi)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);
	
	e = sin(x).*cos(y);
end

function e = dndExtExact(FocalDist,eta,phi)
    x  = FocalDist*cosh(eta).*cos(phi);
    y  = FocalDist*sinh(eta).*sin(phi);
    xn = FocalDist*sinh(eta).*cos(phi);
    yn = FocalDist*cosh(eta).*sin(phi);
    
    e = cos(x).*cos(y).*xn - sin(x).*sin(y).*yn;
end

% function dnde = dndExact(FocalDist,eta,phi,Eta0)
% 	x  = FocalDist*cosh(eta).*cos(phi);
% 	y  = FocalDist*sinh(eta).*sin(phi);
% 	xn = FocalDist*sinh(eta).*cos(phi);
% 	yn = FocalDist*cosh(eta).*sin(phi);
% 	
% 	e = 2*x.*xn - 2.*y.*yn;
% 	
% 	if length(eta)==1 && eta > Eta0
% 		e = cos(x).*cos(y).*xn -  sin(x).*sin(y).*yn;
% 	else
% 		e(eta > Eta0) = cos(x(eta > Eta0)).*cos(y(eta > Eta0)).*xn(eta > Eta0) -  sin(x(eta > Eta0)).*sin(y(eta > Eta0)).*yn(eta > Eta0);
% 	end
% 	
% 	
% end
