function RunLaplacian351

    a=1;%2.5;
    b=1/2;	
	
	x1=-1.1;xn=1.1;
	y1=-0.7;yn=0.7;
	%y1=-1.3;yn=1.3;
	Lx=xn-x1;Ly=yn-y1;
	%ebinf=[];etinf=[];
    ErrIntPre = 0; 	ErrExtPre = 0;	ErrTotPre = 0;
    

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);

    BIn = 1000;
    BOut = 1;

	BType		= 'Fourier';
    
    LinearSolverType = 1;
    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end
    
    
    f   =@(phi) ExtExact(FocalDistance,Eta0,phi);
	g   =@(phi) IntExact(FocalDistance,Eta0,phi);
	gn   =@(phi) dndIntExact(FocalDistance,Eta0,phi);
	
    

	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,g);%(@sin,@sin,1);%(f,dfdn);
	end

    fprintf('Problem 3.51, M=%d, LinearSolverType = %d, BIn=%d, BOut=%d\n', Basis.M, LinearSolverType, BIn, BOut);
    
    
	for n=1:4 %run different grids
		tic
		%build grid
		p=3;%3;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
		
		Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
		ScattererHandle  = @Tools.Scatterer.EllipticScatterer;
		ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',33);

		%------------------------------------------------------------------
		InteriorCoeffsHandle = @Tools.Coeffs.ConstLapCoeffs;
		InteriorCoeffsParams = struct('a',BIn,'b',BIn,'sigma',0);
		
		DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
		%DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
		DiffOpParamsInt = struct('BC_x1',0,'BC_xn', 0,'BC_y1',0,'BC_yn',0, 'LinearSolverType', LinearSolverType);
		%x.^2 - y.^2
		%DiffOpParams = struct(	'BC_x1',Grid.x(1).^2 - Grid.y(2:end-1).^2 ,'BC_xn', Grid.x(end).^2 - Grid.y(2:end-1).^2,...
		%						'BC_y1',(Grid.x(2:end-1).^2 - Grid.y(1).^2).','BC_yn',(Grid.x(2:end-1).^2 - Grid.y(end).^2).' );
		
		IntPrb = Solvers.InteriorHomoLaplacianSolver ...
			(Basis,Grid,InteriorCoeffsHandle,InteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,DiffOp,DiffOpParamsInt);
		
		%------------------------------------------------------------------
		ExteriorCoeffsHandle = @Tools.Coeffs.ConstLapCoeffs;
		ExteriorCoeffsParams = struct('a',BOut,'b',BOut,'sigma',0);
		ExteriorSource          = @Tools.Source.LaplaceSource_IIM351_Exterior;
		SourceParams			= [];
				
		% 		 sin(Grid.x(2:end-1)).*cos(Grid.y(end))
		DiffOpParamsExt = struct('BC_y1', sin(Grid.x(1)).*cos(Grid.y(2:end-1)).','BC_yn',  sin(Grid.x(end)).*cos(Grid.y(2:end-1)).',...
							  'BC_x1',(sin(Grid.x(2:end-1)).*cos(Grid.y(1))),'BC_xn',(sin(Grid.x(2:end-1)).*cos(Grid.y(end))), 'LinearSolverType', LinearSolverType);
		
		ExtPrb =  Solvers.ExteriorLaplacianSolver ...
			(Basis,Grid,ExteriorCoeffsHandle,ExteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,ExteriorSource,SourceParams,DiffOp,DiffOpParamsExt);
							
		%------------------------------------------------------------------
	if 1
		Zeros=spalloc(numel(IntPrb.GridGamma),Basis.NBss,0);
		Eye = speye(Basis.NBss);
		Zeros2=spalloc(Basis.NBss,Basis.NBss,0);
		
		nGG = numel(IntPrb.GridGamma);
		rhs = zeros(2*nGG + 2*Basis.NBss,1);
		%rhs(1:nGG)	= (-IntPrb.TrGF -IntPrb.Qf);
		rhs(nGG+1:2*nGG)= (-ExtPrb.TrGF -ExtPrb.Qf);

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
		
		
		ICu = (x0.^2 - y0.^2) - sin(x0).* cos(y0);
		ICuc = Tools.Basis.FourierBasis.FftCoefs(ICu, Basis.NBss);
		
		ICun = 2*BIn*(x0.*xn0-y0.*yn0) - BOut*(cos(x0).* cos(y0).*xn0 - sin(x0).* sin(y0).*yn0);
		%2*BOut* sin(x0).* cos(y0);
		ICunc = Tools.Basis.FourierBasis.FftCoefs(ICun, Basis.NBss);
		
		rhs(2*nGG+1:2*nGG+Basis.NBss)= ICuc;
		rhs(2*nGG+Basis.NBss+1:end)= ICunc;
		
 		cn = [IntPrb.Q0,IntPrb.Q1,Zeros,Zeros;   
			  Zeros,Zeros,ExtPrb.Q0,ExtPrb.Q1; 
			  Eye,      Zeros2,   -Eye, Zeros2; 
			  Zeros2, BIn*Eye, Zeros2, -BOut*Eye]\rhs;

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
    
		xiex  = IntExact(FocalDistance,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi);
		Intxi(ExtPrb.GridGamma) = xiex; %debig
		
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
    
	
% 	Ex=Tools.Exact.ExLapCrclVarCoeffs(Grid, ExParams);
% 	
%     exact = Ex.u;
	%%%%%%%delete next row
     
	 Intexact = zeros(size(Grid.R));
	 Extexact = zeros(size(Grid.R));
		
	
	Intexact(IntPrb.Scatterer.Np) = IntExact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np));
	Extexact(ExtPrb.Scatterer.Nm) = ExtExact(FocalDistance,ExtPrb.Scatterer.Eta(ExtPrb.Scatterer.Nm),ExtPrb.Scatterer.Phi(ExtPrb.Scatterer.Nm));
	
    ErrInt =norm(Intexact(IntPrb.Scatterer.Np)-Intu(IntPrb.Scatterer.Np),inf);
	%Extetinf(n) =norm(Extexact(ExtPrb.Scatterer.Nm)-Extu(ExtPrb.Scatterer.Nm),inf);
	tmp = Extexact(2:end-1,2:end-1)-Extu(2:end-1,2:end-1);
	ErrExt =norm(tmp(:),inf);
	
    %fprintf('b=%-5.2d,C=%-5.2d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',ExParams.B,ExParams.C,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
    %fprintf('coeffs=%d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',0,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
     
    ErrTot = max(ErrInt,ErrExt);
    
    fprintf('N=%-10dx%-10d\t ErrInt=%d\t rate=%-5.2f ErrExt=%d\t rate=%-5.2f,\t ErrTot=%d\t rate=%-5.2f timeA=%d\n',...
                    Nx,Ny,ErrInt,log2(ErrIntPre/ErrInt),ErrExt,log2(ErrExtPre/ErrExt),ErrTot,log2(ErrTotPre/ErrTot),t1);
    
    ErrIntPre = ErrInt;
    ErrExtPre = ErrExt;
    ErrTotPre = ErrTot;
    
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
