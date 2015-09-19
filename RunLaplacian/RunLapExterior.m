function RunLapExterior
	x1=-1;xn=1;
	y1=-1;yn=1;
	Lx=xn-x1;Ly=yn-y1;
	ebinf=[];etinf=[];
	


	ExParams.B  = 10^3;
	ExParams.C  = 0.1;
	ExParams.r0 = 1/2;
%ExParams.r  = ExParams.r0;

	
	BType		= 'Fourier';

	LinearSolverType = 0;
	if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end
	
    f   =@(theta) Exact(ExParams,theta);
    dfdn=@(theta) drExact(ExParams,theta);
	
	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,dfdn);
	end
	
	for n=1:5 %run different grids
		tic
		%build grid
		p=3;%3;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
		Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
								
		Boundaries.R = [	sqrt(Grid.x(1)^2+ Grid.y(2:end-1).^2).' , sqrt(Grid.x(end)^2+ Grid.y(2:end-1).^2).',  ...
							sqrt(Grid.x(2:end-1).^2+ Grid.y(1).^2).', sqrt(Grid.x(2:end-1).^2+ Grid.y(end).^2).' ];
		
		Ex	= Tools.Exact.ExLapCrclVarCoeffs346(Boundaries, ExParams);		
		DiffOpParams = struct('BC_x1', Ex.u(:,1).','BC_xn', Ex.u(:,2).','BC_y1',Ex.u(:,3),'BC_yn',Ex.u(:,4), 'LinearSolverType', LinearSolverType);
		
		ExtPrb =  Solvers.ExteriorLaplacianSolver( struct(...
                  'Basis'             , Basis, ...
                  'Grid'              , Grid, ...
                  'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                  'CoeffsParams'      , struct('a',ExParams.B,'b',ExParams.B,'sigma',0), ...
                  'ScattererHandle'   , @Tools.Scatterer.PolarScatterer, ...
                  'ScattererParams'   , struct('r0',ExParams.r0,'ExpansionType',33, 'Stencil', 5), ...
                  'CollectRhs'        , CollectRhs, ...
                  'DiffOp'            , @Tools.DifferentialOps.LaplacianOpBCinRhs, ... % @Tools.DifferentialOps.LaplacianOpBCinMat
                  'DiffOpParams'      , DiffOpParams, ...
                  'SourceHandle'      , @Tools.Source.LaplaceSource_IIM346_Exterior, ...
                  'SourceParams'      , ExParams, ...
                  'Extension', @Tools.Extensions.FirstExtension, ...
                  'ExtensionParams',[] ...
                  ));
		
		Q0 = ExtPrb.Q0;
		Q1 = ExtPrb.Q1;
		Qf = ExtPrb.Qf;
		
		TrGF = ExtPrb.TrGF;
		
		cn1 =( Q1 \ ( -Q0*Basis.cn0 - TrGF - Qf)) ;
		
		
		xi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		xi(ExtPrb.GridGamma) = ...
			ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1 + ExtPrb.Wf(ExtPrb.GridGamma);
		
		XiGammaExParams = ExParams;
		XiGammaExParams.r0 = ExtPrb.Scatterer.r;
		xiex = Exact(XiGammaExParams,ExtPrb.Scatterer.th);
		ebinf(n) =norm(xiex -xi(ExtPrb.GridGamma),inf);
		
% 		xi(ExtPrb.GridGamma) = xiex; %test
		u = ExtPrb.P_Omega(xi); 

    %%%%%%%%%%%%%%

		t1=toc;
		
		% % % % % % % % % % % % % % % %
		% Comparison
		% % % % % % % % % % % % % % % %
		
		exact = zeros(Nx,Ny);

		CmprExParams = ExParams;
		CmprExParams.r0 = ExtPrb.Scatterer.R(ExtPrb.Scatterer.Nm);
		exact(ExtPrb.Scatterer.Nm) = Exact(CmprExParams,ExtPrb.Scatterer.Th(ExtPrb.Scatterer.Nm));
		
		t2=toc;
		
		tmp = exact(2:end-1,2:end-1)-u(2:end-1,2:end-1);
		%tmp = exact-u;
		etinf(n) =norm(tmp(:),inf);
		fprintf('b=%-5.2d,C=%-5.2d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',ExParams.B,ExParams.C,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
		
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

