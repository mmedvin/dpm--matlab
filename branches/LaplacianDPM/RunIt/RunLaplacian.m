function RunLaplacian
x1=-1;xn=1;
y1=-1;yn=1;
Lx=xn-x1;Ly=yn-y1;
ebinf=[];etinf=[];



ExParams.B  = 10^-3;
ExParams.C  = 0.1;
ExParams.r0 = 1/2;
%ExParams.r  = ExParams.r0;

    
	BType		= 'Fourier';
    
    f   =@(theta) ExtExact(ExParams,theta);
	g   =@(theta) IntExact(ExParams,theta);
    %dfdn=@(theta) drExact(ExParams,theta);

	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,g,5);%(@sin,@sin,1);%(f,dfdn);
	end

	for n=1:4 %run different grids
		tic
		%build grid
		p=3;%3;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
		
		Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
		SourceParams		= ExParams;
		ScattererClsHandle  = @Tools.Scatterer.PolarScatterer;
		ScattererAddParams  = struct('r0',ExParams.r0,'ExpansionType',33);

		%------------------------------------------------------------------
		InteriorCoeffsClsHandle = @Tools.Coeffs.LaplaceCoeffsPolar;
		InteriorCoeffsClsAddParams = [];
		InteriorSource              = @Tools.Source.LaplaceSourceInterior;
		
		IntPrb = Solvers.InteriorLaplacianSolver ...
			(Basis,Grid,InteriorCoeffsClsHandle,InteriorCoeffsClsAddParams,ScattererClsHandle,ScattererAddParams,InteriorSource,SourceParams);
		
		%------------------------------------------------------------------
		ExteriorCoeffsClsHandle = @Tools.Coeffs.ConstLapCoeffs;
		ExteriorCoeffsAddParams = struct('a',ExParams.B,'b',ExParams.B,'sigma',0);
		ExteriorSource              = @Tools.Source.LaplaceSource;
				
		ExtPrb =  Solvers.ExteriorLaplacianSolver ...
			(Basis,Grid,ExteriorCoeffsClsHandle,ExteriorCoeffsAddParams,ScattererClsHandle,ScattererAddParams,ExteriorSource,SourceParams);
							
		%------------------------------------------------------------------
	if 1
		IntQ = IntPrb.Q;
		ExtQ = ExtPrb.Q;		
	
		rhs = zeros(numel(IntPrb.GridGamma) + numel(ExtPrb.GridGamma),1);
		rhs(1:numel(IntPrb.GridGamma),1)	= (-IntPrb.TrGF -IntPrb.Qf);
		rhs(numel(IntPrb.GridGamma)+1:end,1)= (-ExtPrb.TrGF -ExtPrb.Qf);
	
		cn = [ IntQ ; ExtQ ] \ rhs;
		
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*cn + IntPrb.Wf(IntPrb.GridGamma);
		Intu = IntPrb.P_Omega(Intxi);
				
		Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		Extxi(ExtPrb.GridGamma) = (ExtPrb.W(ExtPrb.GridGamma,:)*cn  + ExtPrb.Wf(ExtPrb.GridGamma));
		Extu = ExtPrb.P_Omega(Extxi);
	else
		%exterior
		Extcn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*Basis.cn0 - ExtPrb.TrGF - ExtPrb.Qf)) ;
		
		
		Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		Extxi(ExtPrb.GridGamma) = ...
			ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*Extcn1 + ExtPrb.Wf(ExtPrb.GridGamma);
		
		%XiGammaExParams = ExParams;
		%XiGammaExParams.r0 = ExtPrb.Scatterer.r;
		%xiex = Exact(XiGammaExParams,ExtPrb.Scatterer.th);
		%ebinf(n) =norm(xiex -xi(ExtPrb.GridGamma),inf);
		
		% % 		xi(ExtPrb.GridGamma) = xiex; %test
		Extu = ExtPrb.P_Omega(Extxi); 
		
		%interior

		Intcn1 =( IntPrb.Q1 \ ( -IntPrb.Q0*Basis.cn0 - IntPrb.TrGF - IntPrb.Qf)) ;
    
        
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = ...
			IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*Intcn1 + IntPrb.Wf(IntPrb.GridGamma);
    
		%XiGammaExParams = ExParams;
		%XiGammaExParams.r0 = IntPrb.Scatterer.r;
		%xiex = Exact(XiGammaExParams,IntPrb.Scatterer.th);
		%ebinf(n) =norm(xiex -xi(IntPrb.GridGamma),inf);

 		%%xi(IntPrb.GridGamma) = xiex; %test
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
     IntCmprExParams = ExParams;
     IntCmprExParams.r0 = Grid.R(IntPrb.Scatterer.Np);
	 
	 Extexact = zeros(size(Grid.R));
	 ExtCmprExParams = ExParams;
     ExtCmprExParams.r0 = Grid.R(ExtPrb.Scatterer.Nm);
		
	Extexact(ExtPrb.Scatterer.Nm) = ExtExact(ExtCmprExParams,ExtPrb.Scatterer.Th(ExtPrb.Scatterer.Nm));
	Intexact(IntPrb.Scatterer.Np) = IntExact(IntCmprExParams,IntPrb.Scatterer.Th(IntPrb.Scatterer.Np));

    t2=toc;
    
    Intetinf(n) =norm(Intexact(IntPrb.Scatterer.Np)-Intu(IntPrb.Scatterer.Np),inf);
	Extetinf(n) =norm(Extexact(ExtPrb.Scatterer.Nm)-Extu(ExtPrb.Scatterer.Nm),inf);
	
    %fprintf('b=%-5.2d,C=%-5.2d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',ExParams.B,ExParams.C,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
    %fprintf('coeffs=%d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',0,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
     
	fprintf('b=%-5.2d,C=%-5.2d,M=%d,N=%-10dx%-10d\t Intetinf=%d\t Extetinf=%d\t timeA=%d\ttimeE=%d\n',ExParams.B,ExParams.C,Basis.M, Nx,Ny,full(Intetinf(n)),full(Extetinf(n)),t1,t2-t1);
end

IntLinf=log2(Intetinf(1:end-1)./Intetinf(2:end))
ExtLinf=log2(Extetinf(1:end-1)./Extetinf(2:end))
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
