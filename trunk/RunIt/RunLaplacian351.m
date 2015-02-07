                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          function RunLaplacian351

    a=1;%2.5;
    b=1/2;	
	
	x1=-2;xn=2;
	%y1=-0.7;yn=0.7;
	%y1=-1.1;yn=1.1;
    
    y1i=-2;yni=2;
	y1e=-2;yne=2;
    
	%Lx=xn-x1;Ly=yn-y1;
	%ebinf=[];etinf=[];
    

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);

    BIn = 1e3;
    BOut = 1e0;

	BType		= 'Fourier';
    
     Order=2; 
     if Order==2 , ExpansionType=33; Stencil=5; else ExpansionType=35; Stencil=9; end
    
    f   =@(phi) ExtExact(FocalDistance,Eta0,phi);
	g   =@(phi) IntExact(FocalDistance,Eta0,phi);
	gn   =@(phi) dndIntExact(FocalDistance,Eta0,phi);
	
    

	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,g);%(@sin,@sin,1);%(f,dfdn);
	end

for	   LinearSolverType = 0
    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end

    % %         ErrIntPre = 0; 	ErrExtPre = 0;	%ErrTotPre = 0; ErrUPre = 0; ErrUrPre = 0; ErrUrrPre = 0;
    % %         ErrUPre = 0; ErrUxPre = 0; ErrUxxPre = 0;  ErrUyPre = 0; ErrUyyPre = 0;  ErrUxyPre = 0;
    
    ErrUInfPre = 0; ErrU2Pre = 0; ErrUxInfPre = 0; ErrUx2Pre = 0; ErrUyInfPre = 0; ErrUy2Pre = 0; ErrUxxInfPre = 0; ErrUxx2Pre = 0; ErrUyyInfPre = 0; ErrUyy2Pre = 0; ErrUxyInfPre = 0; ErrUxy2Pre = 0;
	
    fprintf('Problem 3.51, M=%d, Order = %d,  LinearSolverType = %d, BIn=%d, BOut=%d\n', Basis.M, Order, LinearSolverType, BIn, BOut);
    
    
	for n=1:4 %run different grids
		tic
		%build grid
		p=4;%3;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
		
		%Grid             = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
        
        GridExt                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1e,yne,Ny);
        GridInt                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1i,yni,Ny);
        
		ScattererHandle  = @Tools.Scatterer.EllipticScatterer;
		ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',ExpansionType, 'Stencil', Stencil);

		%------------------------------------------------------------------
		InteriorCoeffsHandle = @Tools.Coeffs.ConstLapCoeffs;
		InteriorCoeffsParams = struct('a',BIn,'b',BIn,'sigma',0);
		
        if Order==4		
            DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat4OrdrConst;%LaplacianOpBCinMat4OrdrConst; %LaplacianOpBCinMat;
        elseif Order==2
            DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
            %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
        end
		DiffOpParamsInt = struct('BC_x1',0,'BC_xn', 0,'BC_y1',0,'BC_yn',0, 'LinearSolverType', LinearSolverType);
		%x.^2 - y.^2
		%DiffOpParams = struct(	'BC_x1',Grid.x(1).^2 - Grid.y(2:end-1).^2 ,'BC_xn', Grid.x(end).^2 - Grid.y(2:end-1).^2,...
		%						'BC_y1',(Grid.x(2:end-1).^2 - Grid.y(1).^2).','BC_yn',(Grid.x(2:end-1).^2 - Grid.y(end).^2).' );
		
		IntPrb = Solvers.InteriorHomoLaplacianSolver ...
			(Basis,GridInt,InteriorCoeffsHandle,InteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,DiffOp,DiffOpParamsInt);
		
		%------------------------------------------------------------------
		ExteriorCoeffsHandle = @Tools.Coeffs.ConstLapCoeffs;
		ExteriorCoeffsParams = struct('a',BOut,'b',BOut,'sigma',0);
		ExteriorSource          = @Tools.Source.LaplaceSource_IIM351_Exterior;
		SourceParams			= [];
				
		% 		 sin(Grid.x(2:end-1)).*cos(Grid.y(end))
        if isequal(DiffOp , @Tools.DifferentialOps.LaplacianOpBCinRhs)
            DiffOpParamsExt = struct('BC_y1', sin(GridExt.x(1)).*cos(GridExt.y(2:end-1)).','BC_yn',  sin(GridExt.x(end)).*cos(GridExt.y(2:end-1)).',...
                'BC_x1',(sin(GridExt.x(2:end-1)).*cos(GridExt.y(1))),'BC_xn',(sin(GridExt.x(2:end-1)).*cos(GridExt.y(end))), 'LinearSolverType', LinearSolverType);
        elseif isequal(DiffOp , @Tools.DifferentialOps.LaplacianOpBCinMat) || isequal(DiffOp , @Tools.DifferentialOps.LaplacianOpBCinMat4OrdrConst)
            DiffOpParamsExt = struct('BC_y1', sin(GridExt.x(1)).*cos(GridExt.y).','BC_yn',  sin(GridExt.x(end)).*cos(GridExt.y).',...
                'BC_x1',(sin(GridExt.x).*cos(GridExt.y(1))),'BC_xn',(sin(GridExt.x).*cos(GridExt.y(end))), 'LinearSolverType', LinearSolverType);
        end
        
		
		ExtPrb =  Solvers.ExteriorLaplacianSolver... ExteriorLaplacianSolver ... ExtLapSolver4Ordr
			(Basis,GridExt,ExteriorCoeffsHandle,ExteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,ExteriorSource,SourceParams,DiffOp,DiffOpParamsExt);
							
		%------------------------------------------------------------------
	if 1
        nGGext = numel(ExtPrb.GridGamma);
        nGGint = numel(IntPrb.GridGamma);
        
		ZerosInt=spalloc(nGGint,Basis.NBss,0);
        ZerosExt=spalloc(nGGext,Basis.NBss,0);
        
		Eye = speye(Basis.NBss);
		Zeros2=spalloc(Basis.NBss,Basis.NBss,0);
		
		rhs = zeros(nGGext+nGGint + 2*Basis.NBss,1);
		%rhs(1:nGG)	= (-IntPrb.TrGF -IntPrb.Qf);
		rhs((nGGint+1):(nGGext+nGGint))= -(ExtPrb.TrGF +ExtPrb.Qf);

		% [u_in] = x0^2 - y0^2 
		% [u_out] = sin(x0).* cos(y0) 
		% [ beta un_in] = 0
		% [beta un_out] = - 2 b sin(x0).* cos(y0)
		
		phi=linspace(0,2*pi, Basis.NBss+1);
		phi = phi(1:Basis.NBss);

		x0 = FocalDistance*cosh(Eta0).*cos(phi);
		y0 = FocalDistance*sinh(Eta0).*sin(phi);
		xn0 = FocalDistance*sinh(Eta0).*cos(phi);
		yn0 = FocalDistance*cosh(Eta0).*sin(phi);
		
		if 0
		ICu = (x0.^2 - y0.^2) - sin(x0).* cos(y0);
		ICuc = Tools.Basis.FourierBasis.FftCoefs(ICu, Basis.NBss);
		
		ICun = 2*BIn*(x0.*xn0-y0.*yn0) - BOut*(cos(x0).* cos(y0).*xn0 - sin(x0).* sin(y0).*yn0);
		%2*BOut* sin(x0).* cos(y0);
		ICunc = Tools.Basis.FourierBasis.FftCoefs(ICun, Basis.NBss);
		
		rhs((nGGext+nGGint+1):(nGGext+nGGint+Basis.NBss))= ICuc;
		rhs((nGGext+nGGint+Basis.NBss+1):end)= ICunc;
		
 		cn = [IntPrb.Q0,IntPrb.Q1,ZerosInt,ZerosInt;   
			  ZerosExt,ZerosExt,ExtPrb.Q0,ExtPrb.Q1; 
			  Eye,      Zeros2,   -Eye, Zeros2; 
			  Zeros2, BIn*Eye, Zeros2, -BOut*Eye]\rhs;
        else
            sinxpy = sin(x0+y0);
            sinxmy = sin(x0-y0);
            cosxpy = cos(x0+y0);
            cosxmy = cos(x0-y0);
            cosxcosy = (cosxmy + cosxpy)/2;
            sinxsiny = (cosxmy - cosxpy)/2;
            sinxcosy = (sinxpy + sinxmy)/2;
            
            ICuin = (x0.^2 - y0.^2) ;
            ICuout= sinxcosy; %sin(x0).* cos(y0);
            ICuinc = Tools.Basis.FourierBasis.FftCoefs(ICuin, Basis.NBss);
            ICuoutc = Tools.Basis.FourierBasis.FftCoefs(ICuout, Basis.NBss);
            
            %ICuinn = 2*BIn*(x0.*xn0-y0.*yn0) ;
            %ICuoutn = BOut*(cos(x0).* cos(y0).*xn0 - sin(x0).* sin(y0).*yn0);
            ICuinn = 2*(x0.*xn0-y0.*yn0) ;
            ICuoutn = (cosxcosy.*xn0 - sinxsiny.*yn0);%(cos(x0).* cos(y0).*xn0 - sin(x0).* sin(y0).*yn0);
            %2*BOut* sin(x0).* cos(y0);
            ICuinnc = Tools.Basis.FourierBasis.FftCoefs(ICuinn, Basis.NBss);
            ICuoutnc = Tools.Basis.FourierBasis.FftCoefs(ICuoutn, Basis.NBss);

            rhs((nGGext+nGGint+1):(nGGext+nGGint+Basis.NBss))= ICuinc(:);
            rhs((nGGext+nGGint+Basis.NBss+1):(nGGext+nGGint+2*Basis.NBss))= ICuoutc;
            rhs((nGGext+nGGint+2*Basis.NBss+1):(nGGext+nGGint+3*Basis.NBss))= ICuinnc;
            rhs((nGGext+nGGint+3*Basis.NBss+1):(nGGext+nGGint+4*Basis.NBss))= ICuoutnc;
            
            A = [IntPrb.Q0,IntPrb.Q1,ZerosInt,ZerosInt;
                ZerosExt,ZerosExt,ExtPrb.Q0,ExtPrb.Q1;
                Eye,      Zeros2, Zeros2  , Zeros2;
                Zeros2,      Zeros2,   Eye, Zeros2;
                Zeros2, Eye, Zeros2, Zeros2;
                Zeros2, Zeros2, Zeros2, Eye];
               % Zeros2, BIn*Eye, Zeros2, Zeros2;
               % Zeros2, Zeros2, Zeros2, BOut*Eye];
           % [cn,flag] = lsqr(A'*A,A'*rhs,1e-14,Nx);
            %cn = (A'*A)\(A'*rhs);
            
           % [Q,R,P]=qr(A'*A);
            %cn = P*( R\(Q\(A'*rhs)));
            cn=A\rhs;
        end
		Intcn=cn(1:2*Basis.NBss);
		Extcn= cn(2*Basis.NBss+1:end); 
		
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*Intcn;% + IntPrb.Wf(IntPrb.GridGamma);
		Intu = IntPrb.P_Omega(Intxi);
				
		Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		Extxi(ExtPrb.GridGamma) = (ExtPrb.W(ExtPrb.GridGamma,:)*Extcn  + ExtPrb.Wf(ExtPrb.GridGamma));
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
			ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*Extcn1 + ExtPrb.Wf(ExtPrb.GridGamma);
		
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
    
    
    [ErrUInf,ErrU2] = cmpr([Intu(IntPrb.Scatterer.Mp);Extu(ExtPrb.Scatterer.Mm)],[Intexact(IntPrb.Scatterer.Mp);Extexact(ExtPrb.Scatterer.Mm)]);
    
    CDExt = Tools.Common.SecondDerivative(GridExt.Nx,GridExt.Ny,GridExt.dx,GridExt.dy);
    CDInt = Tools.Common.SecondDerivative(GridInt.Nx,GridInt.Ny,GridInt.dx,GridInt.dy);

    [eTux,eTuy,eTuxx,eTuyy,eTuxy]       = CDExt.CartesianDerivatives(Extu);
    [iTux,iTuy,iTuxx,iTuyy,iTuxy]       = CDInt.CartesianDerivatives(Intu);
    [eExux,eExuy,eExuxx,eExuyy,eExuxy]  = CDExt.CartesianDerivatives(Extexact);
    [iExux,iExuy,iExuxx,iExuyy,iExuxy]  = CDInt.CartesianDerivatives(Intexact);
    
    [ErrUxInf,ErrUx2]   = cmpr([iTux(IntPrb.Scatterer.Mp);eTux(ExtPrb.Scatterer.Mm)],[iExux(IntPrb.Scatterer.Mp);eExux(ExtPrb.Scatterer.Mm)]);
    [ErrUyInf,ErrUy2]   = cmpr([iTuy(IntPrb.Scatterer.Mp);eTuy(ExtPrb.Scatterer.Mm)],[iExuy(IntPrb.Scatterer.Mp);eExuy(ExtPrb.Scatterer.Mm)]);
	
    [ErrUxxInf,ErrUxx2] = cmpr([iTuxx(IntPrb.Scatterer.Mp);eTuxx(ExtPrb.Scatterer.Mm)],[iExuxx(IntPrb.Scatterer.Mp);eExuxx(ExtPrb.Scatterer.Mm)]);
    [ErrUyyInf,ErrUyy2] = cmpr([iTuyy(IntPrb.Scatterer.Mp);eTuyy(ExtPrb.Scatterer.Mm)],[iExuyy(IntPrb.Scatterer.Mp);eExuyy(ExtPrb.Scatterer.Mm)]);
    [ErrUxyInf,ErrUxy2] = cmpr([iTuxy(IntPrb.Scatterer.Mp);eTuxy(ExtPrb.Scatterer.Mm)],[iExuxy(IntPrb.Scatterer.Mp);eExuxy(ExtPrb.Scatterer.Mm)]);
%------------------------------------------------------------------
     
    fprintf('N=%-6dx%-7d Eu=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Eux=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euxx=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euyy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euxy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f timeA=%-6.2f\n',...
        Nx,Ny,ErrUInf,   ErrU2,   log2(ErrUInfPre/ErrUInf),     log2(ErrU2Pre/ErrU2), ...
              ErrUxInf,  ErrUx2,  log2(ErrUxInfPre/ErrUxInf),   log2(ErrUx2Pre/ErrUx2), ...  
              ErrUyInf,  ErrUy2,  log2(ErrUyInfPre/ErrUyInf),   log2(ErrUy2Pre/ErrUy2), ...  
              ErrUxxInf, ErrUxx2, log2(ErrUxxInfPre/ErrUxxInf), log2(ErrUxx2Pre/ErrUxx2), ...
              ErrUyyInf, ErrUyy2, log2(ErrUyyInfPre/ErrUyyInf), log2(ErrUyy2Pre/ErrUyy2), ...  
              ErrUxyInf, ErrUxy2, log2(ErrUxyInfPre/ErrUxyInf), log2(ErrUxy2Pre/ErrUxy2), ...
            t1);
		
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



function [Linf,L2] = cmpr(ex,u,GG)

    if nargin==2, GG=[];end
    tmp = ex - u;
    
    if ~isempty(GG)
        tmp(GG) = 0;
        u(GG)=0;
    end

    Linf = norm(tmp(:),inf)/norm(u(:),inf);
    L2   = norm(tmp(:),2)/norm(u(:),2);
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
