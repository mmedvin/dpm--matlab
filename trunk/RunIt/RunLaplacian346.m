function RunLaplacian346
x1=-1;xn=1;
y1=-1;yn=1;
Lx=xn-x1;Ly=yn-y1;
%ebinf=[];etinf=[];
	


ExParams.B  = 10^(-3);
ExParams.C  = 0.1;
ExParams.r0 = 1/2;

    
	BType		= 'Fourier';
	    
    f   =@(theta) ExtExact(ExParams,theta);
	g   =@(theta) IntExact(ExParams,theta);
    %dfdn=@(theta) drExact(ExParams,theta);

	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,g);%(@sin,@sin,1);%(f,dfdn);
	end

	
for	LinearSolverType = 0
	if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end

	ErrIntPre = 0; 	ErrExtPre = 0;	ErrTotPre = 0; ErrUPre = 0; ErrUrPre = 0; ErrUrrPre = 0;
	
	fprintf('Problem 3.46, M=%d, LinearSolverType = %d, r0=%-4.2f, C=%-6.3d,B=%-6.3d\n', Basis.M, LinearSolverType, ExParams.r0,ExParams.C,ExParams.B);
	
	for n=1:5 %run different grids
		tic
		%build grid
		p=3;%3;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
		
		Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
		SourceParams		= ExParams;
		ScattererHandle  = @Tools.Scatterer.PolarScatterer;
		ScattererParams  = struct('r0',ExParams.r0,'ExpansionType',33);

		%------------------------------------------------------------------
		InteriorCoeffsHandle = @Tools.Coeffs.LaplaceCoeffsPolar;
		InteriorCoeffsParams = [];
		InteriorSource              = @Tools.Source.LaplaceSource_IIM346_Interior;
		
		DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
		
		%Boundaries.R = [	sqrt(Grid.x(1)^2+ Grid.y(2:end-1).^2).' , sqrt(Grid.x(end)^2+ Grid.y(2:end-1).^2).',  ...
		%					sqrt(Grid.x(2:end-1).^2+ Grid.y(1).^2).', sqrt(Grid.x(2:end-1).^2+ Grid.y(end).^2).' ];
		Boundaries.R = [	sqrt(Grid.x(2:end-1).^2+ Grid.y(1).^2).', sqrt(Grid.x(2:end-1).^2+ Grid.y(end).^2).',  ...
							sqrt(Grid.x(1)^2+ Grid.y(2:end-1).^2).' , sqrt(Grid.x(end)^2+ Grid.y(2:end-1).^2).'	];
		
		
		
		Exact	= Tools.Exact.ExLapCrclVarCoeffs346(Boundaries, SourceParams);		
		DiffOpParams = struct('BC_x1', Exact.u(:,1).','BC_xn', Exact.u(:,2).','BC_y1',Exact.u(:,3),'BC_yn',Exact.u(:,4), 'LinearSolverType', LinearSolverType);
				
		IntPrb = Solvers.InteriorLaplacianSolver ...
			(Basis,Grid,InteriorCoeffsHandle,InteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,InteriorSource,SourceParams,DiffOp,DiffOpParams);
		
		%------------------------------------------------------------------
		ExteriorCoeffsHandle = @Tools.Coeffs.ConstLapCoeffs;
		ExteriorCoeffsParams = struct('a',ExParams.B,'b',ExParams.B,'sigma',0);
		ExteriorSource              = @Tools.Source.LaplaceSource_IIM346_Exterior;
				
		ExtPrb =  Solvers.ExteriorLaplacianSolver ...
			(Basis,Grid,ExteriorCoeffsHandle,ExteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,ExteriorSource,SourceParams,DiffOp,DiffOpParams);
							
		%------------------------------------------------------------------
	if 1
		Zeros=spalloc(numel(IntPrb.GridGamma),Basis.NBss,0);
		Eye = speye(Basis.NBss);
		Zeros2=spalloc(Basis.NBss,Basis.NBss,0);
		
		nGG = numel(IntPrb.GridGamma);
		rhs = zeros(2*nGG + Basis.NBss,1);
		rhs(1:nGG)	= (-IntPrb.TrGF -IntPrb.Qf);
		rhs(nGG+1:2*nGG)= (-ExtPrb.TrGF -ExtPrb.Qf);
		
		IC = - (ExParams.C + 2*(ExParams.r0^2).*(1 - ExParams.B + ExParams.r0.^2))./(ExParams.B.*ExParams.r0);
		ICc = Tools.Basis.FourierBasis.FftCoefs(IC*ones(Basis.NBss,1), Basis.NBss);
		
		rhs(2*nGG+1:end)= ICc;
		
 		cn = [IntPrb.Q0,IntPrb.Q1,Zeros;   ExtPrb.Q0,Zeros,ExtPrb.Q1; Zeros2, Eye, -Eye]\rhs;

		Intcn=cn(1:2*Basis.NBss);
		Extcn=[cn(1:Basis.NBss); cn(2*Basis.NBss+1:end)]; 
		
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*Intcn + IntPrb.Wf(IntPrb.GridGamma);
		Intu = IntPrb.P_Omega(Intxi);
				
		Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		Extxi(ExtPrb.GridGamma) = (ExtPrb.W(ExtPrb.GridGamma,:)*Extcn  + ExtPrb.Wf(ExtPrb.GridGamma));

		Extu = spalloc(Nx,Ny,numel(ExtPrb.Scatterer.Nm));
		tmp = ExtPrb.P_Omega(Extxi);
		Extu(ExtPrb.Scatterer.Nm) = tmp(ExtPrb.Scatterer.Nm);
		
		
	else
		%exterior
		Extcn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*Basis.cn0 - ExtPrb.TrGF - ExtPrb.Qf)) ;
		
		
		Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		Extxi(ExtPrb.GridGamma) = ...
			ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*Extcn1 + ExtPrb.Wf(ExtPrb.GridGamma);
		
		Extu = ExtPrb.P_Omega(Extxi); 
		
		%interior

		Intcn1 =( IntPrb.Q1 \ ( -IntPrb.Q0*Basis.cn0 - IntPrb.TrGF - IntPrb.Qf)) ;
    
        
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = ...
			IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*Intcn1 + IntPrb.Wf(IntPrb.GridGamma);
    
		Intu = IntPrb.P_Omega(Intxi);
	end
     
   t1=toc; 

%------------------------------------------------------------------
% Comparison
%------------------------------------------------------------------
    
	
 	Ex=Tools.Exact.ExLapCrclVarCoeffs346(Grid, ExParams);

	u=Ex.u; %zeros(size(Grid.R));
	%u(ExtPrb.Scatterer.Mm) = Extu(ExtPrb.Scatterer.Mm);
	u(2:end-1,2:end-1) = Extu(2:end-1,2:end-1);
	u(IntPrb.Scatterer.Mp) = Intu(IntPrb.Scatterer.Mp);
	ErrU = norm(Ex.u(:) - u(:),inf);

	dxdr = cos(Grid.Theta);
	dydr = sin(Grid.Theta);
	
	[ux,uy] = gradient(u,Grid.dx,Grid.dy);
	ur = Ex.dudr;
	tmp = ux.*dxdr + uy.*dydr;
	ur(2:end-1,2:end-1) = tmp(2:end-1,2:end-1); 
	tmp = Ex.dudr - ur;
	tmp(IntPrb.Scatterer.GridGamma) = 0;
	ErrUr = norm(tmp(:),inf);
	
	%[urx,ury] = gradient(ur,Grid.dx,Grid.dy);
	[uxx,uxy] = gradient(ux,Grid.dx,Grid.dy);
	[uyx,uyy] = gradient(uy,Grid.dx,Grid.dy);
	
	tmp=uxx.*dxdr.^2 + 2*uxy.*dxdr.*dydr + uyy.*dydr.^2;
	
	urr = Ex.d2udr2;
	%tmp = urx.*dxdr + ury.*dydr;
	urr(3:end-2,3:end-2) = tmp(3:end-2,3:end-2); 
	tmp = Ex.d2udr2 - urr;
	tmp(IntPrb.Scatterer.GridGamma) = 0;
	%ErrUrr  = norm(tmp(Grid.R<0.78 & Grid.R>0.73),inf);
	%tmp( Grid.R>0.35  )=0;
	tmp( Grid.R<0.7 )=0;
	ErrUrr  = norm(tmp( :),inf);
%------------------------------------------------------------------
	
	
	 Intexact = zeros(size(Grid.R));
     IntCmprExParams = ExParams;
     IntCmprExParams.r0 = Grid.R(IntPrb.Scatterer.Np);
	 
	 Extexact = zeros(size(Grid.R));
	 ExtCmprExParams = ExParams;
     ExtCmprExParams.r0 = Grid.R(ExtPrb.Scatterer.Nm);
		
	Extexact(ExtPrb.Scatterer.Nm) = ExtExact(ExtCmprExParams,ExtPrb.Scatterer.Th(ExtPrb.Scatterer.Nm));
	Intexact(IntPrb.Scatterer.Np) = IntExact(IntCmprExParams,IntPrb.Scatterer.Th(IntPrb.Scatterer.Np));
    
    ErrInt =norm(Intexact(IntPrb.Scatterer.Np)-Intu(IntPrb.Scatterer.Np),inf);
	%Extetinf =norm(Extexact(ExtPrb.Scatterer.Nm)-Extu(ExtPrb.Scatterer.Nm),inf);
	tmp=Extexact(2:end-1,2:end-1)-Extu(2:end-1,2:end-1);
	ErrExt =norm(tmp(:),inf);
	
	ErrTot = max(ErrInt,ErrExt);
	
    %fprintf('b=%-5.2d,C=%-5.2d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',ExParams.B,ExParams.C,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
    %fprintf('coeffs=%d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',0,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
     
	fprintf('N=%-6dx%-7d EInt=%-10.4d rt=%-6.2f EExt=%-10.4d rt=%-6.2f ETot=%d\t rt=%-6.2f EU=%d\t rt=%-6.2f EUr=%d\t rt=%-6.2f EUrr=%d\t rt=%-6.2f  timeA=%-6.2f\n',...
		Nx,Ny,ErrInt,log2(ErrIntPre/ErrInt),ErrExt,log2(ErrExtPre/ErrExt),ErrTot,log2(ErrTotPre/ErrTot),ErrU,log2(ErrUPre/ErrU),ErrUr,log2(ErrUrPre/ErrUr),ErrUrr,log2(ErrUrrPre/ErrUrr),t1);
	
	ErrIntPre = ErrInt;
	ErrExtPre = ErrExt;
	ErrTotPre = ErrTot;
	ErrUPre	  = ErrU;
	ErrUrPre  = ErrUr;
	ErrUrrPre  = ErrUrr;
	
end
end
%IntLinf=log2(Intetinf(1:end-1)./Intetinf(2:end))
%ExtLinf=log2(Extetinf(1:end-1)./Extetinf(2:end))
% Lbinf=log2(ebinf(1:end-1)./ebinf(2:end))

end
	
function e = ExtExact(Params,theta)
    r  = Params.r0   ;
	if size(r)==1
		r = ones(size(theta)).*r;
	end
	e     = (1 - 1/8/Params.B - 1/Params.B)/4 + ( (r.^4)/2 + r.^2 )/Params.B + Params.C*log(2*r)/Params.B;
end
	
function e = IntExact(Params,theta)
	r = Params.r0;
	if size(r)==1
		r = ones(size(theta)).*r;
	end
	
e = r.^2;
end
