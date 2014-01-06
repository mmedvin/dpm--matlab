function RunLapExterior
	x1=-1;xn=1;
	y1=-1;yn=1;
	Lx=xn-x1;Ly=yn-y1;
	ebinf=[];etinf=[];
	
	ExParams.B  = 10^3;
	ExParams.C  = 0.1;
	ExParams.r0 = 1/4;
%ExParams.r  = ExParams.r0;

	
	BType		= 'Fourier';
    
    f   =@(theta) Exact(ExParams,theta);
    dfdn=@(theta) drExact(ExParams,theta);
	
	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,dfdn);
	end
	
	for n=1:4 %run different grids
		tic
		%build grid
		p=4;%3;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
		Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
		
		WaveNumberClsHandle = @Tools.Coeffs.ConstLapCoeffs;%WaveNumberElliptical;%ConstantWaveNumber;%WaveNumberPolarR
		WaveNumberAddParams = struct('a',ExParams.B,'b',ExParams.B,'sigma',0);
		
		ScattererClsHandle  = @Tools.Scatterer.PolarScatterer;
		ScattererAddParams  = struct('r0',ExParams.r0,'ExpansionType',33);
		
		Source              = @Tools.Source.LaplaceSource;
		SourceParams		= ExParams;
		
		IntPrb =  Solvers.InteriorLaplacianSolver ...
			(Basis,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams,Source,SourceParams);
		
		Q0 = IntPrb.Q0;
		Q1 = IntPrb.Q1;
		Qf = IntPrb.Qf;
		
		TrGF = IntPrb.TrGF;
		
		cn1 =( Q1 \ ( -Q0*Basis.cn0 - TrGF - Qf)) ;
		
		
		xi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		xi(IntPrb.GridGamma) = ...
			IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1 + IntPrb.Wf(IntPrb.GridGamma);
		
		ebinf(n)=0;
		
		XiGammaExParams = ExParams;
		XiGammaExParams.r0 = IntPrb.Scatterer.r;
		xiex = Exact(XiGammaExParams,IntPrb.Scatterer.th);
		ebinf(n) =norm(xiex -xi(IntPrb.GridGamma),inf);
		
% 		xi(IntPrb.GridGamma) = xiex; %test
		u = IntPrb.P_Omega(xi);
		
		t1=toc;
		
		% % % % % % % % % % % % % % % %
		% Comparison
		% % % % % % % % % % % % % % % %
		
		
		exact = zeros(Nx,Ny);

		CmprExParams = ExParams;
		CmprExParams.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
		exact(IntPrb.Scatterer.Np) = Exact(CmprExParams,IntPrb.Scatterer.Th(IntPrb.Scatterer.Np));
		
		t2=toc;
		
		etinf(n) =norm(exact(:)-u(:),inf);
		fprintf('b=%-5.2d,C=%-5.2d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',ExParams.B,ExParams.C,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
		
		
		% figure, plot(x,abs(u(fix(Ny/2),:)),'r+-',x,abs(exact(fix(Ny/2),:)),'bo');
		%     saveas(gcf,['inhomo_k=' num2str(k) 'N=' num2str(Nx) '.jpg'],'jpg');
		% plot(x,abs(u(:,fix(Ny/2))),'r+-',x,abs(exact(:,fix(Ny/2))),'bo');
		
	end
	
	Linf=log2(etinf(1:end-1)./etinf(2:end))
	Lbinf=log2(ebinf(1:end-1)./ebinf(2:end))
	
end



function e = Exact(Params,theta)
    r  = Params.r0   ;
	if size(r)==1
		r = ones(size(theta)).*r;
	end
	e     = (1 - 1/8/Params.B - 1/Params.B)/4 + ( (r.^4)/2 + r.^2 )/Params.B + Params.C*log(2*r)/Params.B;
end

function er = drExact(Params,theta)
    r  = Params.r0   ;
	if size(r)==1
		r = ones(size(theta)).*r;
	end
  
	er = 2*( r.^3 + r )/Params.B + Params.C/Params.B./r;
end

