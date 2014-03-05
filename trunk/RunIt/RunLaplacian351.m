function RunLaplacian351

    a=1;%2.5;
    b=1/2;	
	
	x1=-1.2;xn=1.2;
	y1=-0.3;yn=0.3;
	Lx=xn-x1;Ly=yn-y1;
	ebinf=[];etinf=[];

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);

BIn = 1000;
BOut = 1;

	BType		= 'Fourier';
    
    f   =@(phi) ExtExact(FocalDistance,Eta0,phi);
	g   =@(phi) IntExact(FocalDistance,Eta0,phi);
    

	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,g);%(@sin,@sin,1);%(f,dfdn);
	end

	for n=1:5 %run different grids
		tic
		%build grid
		p=3;%3;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
		
		Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
		ScattererClsHandle  = @Tools.Scatterer.EllipticScatterer;
		ScattererAddParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',33);

		%------------------------------------------------------------------
		InteriorCoeffsClsHandle = @Tools.Coeffs.ConstLapCoeffs;
		InteriorCoeffsClsAddParams = struct('a',BIn,'b',BIn,'sigma',0);
		
		IntPrb = Solvers.InteriorHomoLaplacianSolver ...
			(Basis,Grid,InteriorCoeffsClsHandle,InteriorCoeffsClsAddParams,ScattererClsHandle,ScattererAddParams);
		
		%------------------------------------------------------------------
		ExteriorCoeffsClsHandle = @Tools.Coeffs.ConstLapCoeffs;
		ExteriorCoeffsAddParams = struct('a',BOut,'b',BOut,'sigma',0);
		ExteriorSource          = @Tools.Source.LaplaceSource_IIM351_Exterior;
		SourceParams			= [];
				
		ExtPrb =  Solvers.ExteriorLaplacianSolver ...
			(Basis,Grid,ExteriorCoeffsClsHandle,ExteriorCoeffsAddParams,ScattererClsHandle,ScattererAddParams,ExteriorSource,SourceParams);
							
		%------------------------------------------------------------------
	if 1
		
		%x0 = FocalDist*cosh(Eta0).*cos(phi);
		%y0 = FocalDist*sinh(Eta0).*sin(phi);
		% [u_in] = x0^2 - y0^2 
		% [u_out] = sin(x0).* cos(y0) 
		% [ beta un_in] = 0
		% [beta un_out] = - 2 b sin(x0).* cos(y0)
		
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
		Extu = ExtPrb.P_Omega(Extxi);
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
	
function e = IntExact(FocalDist,eta,phi)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);
	
	e = x.^2 - y.^2;	
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
