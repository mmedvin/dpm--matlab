function RunLapBL54
% 5.3  bealy-layton

    %a=0.9; b=0.7; %table 1
    a=0.9; b=0.1; %table 2
    
	x1=-1.1;xn=1.1;
	y1=-1.1;yn=1.1;

	Lx=xn-x1; Ly=yn-y1;

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
    
    ExParams.c=9;
    ExParams.d=8;

	BType		= 'Fourier';
    
	f   =@(phi) Exact(FocalDistance,Eta0,phi,ExParams);
	fn   =@(phi) dndExact(FocalDistance,Eta0,phi,ExParams);
	    
	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,fn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,fn);
	end

for	   LinearSolverType = 0
    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end

    ErrUInfPre = 0; ErrU2Pre = 0; ErrUxInfPre = 0; ErrUx2Pre = 0; ErrUyInfPre = 0; ErrUy2Pre = 0; ErrUxxInfPre = 0; ErrUxx2Pre = 0; ErrUyyInfPre = 0; ErrUyy2Pre = 0; ErrUxyInfPre = 0; ErrUxy2Pre = 0;
	
    fprintf('Problem BL53, M=%d, LinearSolverType = %d \n', Basis.M, LinearSolverType);
        
	for n=1:4 %6 %run different grids
		tic
		%build grid
        p=20;
		Nx=p*2^n; Ny=p*2^n;
		
		Grid             = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
		ScattererHandle  = @Tools.Scatterer.EllipticScatterer;
		ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',33);

		%------------------------------------------------------------------

		CoeffsHandle    = @Tools.Coeffs.ConstLapCoeffs;
		CoeffsParams    = struct('a',1,'b',1,'sigma',0);
		Source          = @Tools.Source.LaplaceSource_BL54_Exterior;
		SourceParams	= ExParams;
		
        DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
		%DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
               
		DiffOpParamsExt = struct(   'BC_y1',  (Grid.x(1).^ExParams.c)       .*(Grid.y(2:end-1).^ExParams.d).', ...
                                    'BC_yn',  (Grid.x(end).^ExParams.c)     .*(Grid.y(2:end-1).^ExParams.d).',...
                                    'BC_x1',( (Grid.x(2:end-1).^ExParams.c) .*(Grid.y(1).^ExParams.d)       ), ...
                                    'BC_xn',( (Grid.x(2:end-1).^ExParams.c) .*(Grid.y(end).^ExParams.d)     ), 'LinearSolverType', LinearSolverType);
		
		ExtPrb =  Solvers.ExteriorLaplacianSolver ...
			(Basis,Grid,CoeffsHandle,CoeffsParams,ScattererHandle,ScattererParams,CollectRhs,Source,SourceParams,DiffOp,DiffOpParamsExt);
							
		%------------------------------------------------------------------

		%exterior
		
		EQ0 = ExtPrb.Q0;
		EQ1 = ExtPrb.Q1;
		ETrGF = ExtPrb.TrGF;
		EQf = ExtPrb.Qf;
		Extcn1 =( EQ1 \ ( -EQ0*Basis.cn0 -ETrGF - EQf)) ;
		
		Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		Extxi(ExtPrb.GridGamma) = ...
			ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*Extcn1 + ExtPrb.Wf(ExtPrb.GridGamma);
		
		%xiex  = ExtExact(FocalDistance,ExtPrb.Scatterer.eta,ExtPrb.Scatterer.phi);
		%Extxi(ExtPrb.GridGamma) = xiex; %debig
		
		Extu = ExtPrb.P_Omega(Extxi); 
		   
        t1=toc;

   %------------------------------------------------------------------
    % Comparison
    %------------------------------------------------------------------
    
    exact = zeros(size(Grid.R));
    exact(ExtPrb.Scatterer.Nm) = Exact(FocalDistance,ExtPrb.Scatterer.Eta(ExtPrb.Scatterer.Nm),ExtPrb.Scatterer.Phi(ExtPrb.Scatterer.Nm),ExParams);    
    
    u = zeros(size(Grid.R));
    u(ExtPrb.Scatterer.Nm) = Extu(ExtPrb.Scatterer.Nm);
    
    u(1,:)   = exact(1,:);
    u(end,:) = exact(end,:);
    u(:,1)   = exact(:,1);
    u(:,end) = exact(:,end);
    
    
    %u = Exact(FocalDistance,ExtPrb.Scatterer.Eta,ExtPrb.Scatterer.Phi);
    %u(ExtPrb.Scatterer.Mm) = Extu(ExtPrb.Scatterer.Mm);
	%u(2:end-1,2:end-1) = Extu(2:end-1,2:end-1);
	%u(ExtPrb.Scatterer.Mp) = 0;
	[ErrUInf,ErrU2] = cmpr(exact,u,0);
    
    %ErrUInf = norm(exact(:) - u(:),inf)/norm(exact(:),inf);
    %ErrU2 = norm(exact(:) - u(:),2)/norm(u(:),2);
 
    CD = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy);
	[Tux,Tuy,Tuxx,Tuyy,Tuxy] = CD.CartesianDerivatives(u(:));
    
    Tux  = reshape(Tux,Grid.Nx,Grid.Ny);
    Tuy  = reshape(Tuy,Grid.Nx,Grid.Ny);
    Tuxx = reshape(Tuxx,Grid.Nx,Grid.Ny);
    Tuyy = reshape(Tuyy,Grid.Nx,Grid.Ny);
    Tuxy = reshape(Tuxy,Grid.Nx,Grid.Ny);
    
    [Exux,Exuy,Exuxx,Exuyy,Exuxy] = CD.CartesianDerivatives(exact(:));
    
    Exux  = reshape(Exux,Grid.Nx,Grid.Ny);
    Exuy  = reshape(Exuy,Grid.Nx,Grid.Ny);
    Exuxx = reshape(Exuxx,Grid.Nx,Grid.Ny);
    Exuyy = reshape(Exuyy,Grid.Nx,Grid.Ny);
    Exuxy = reshape(Exuxy,Grid.Nx,Grid.Ny);
        
    ux = Exux;
    ux(2:end-1,2:end-1) = Tux(2:end-1,2:end-1);
    [ErrUxInf,ErrUx2] = cmpr(Exux,ux,1);
    
    uy = Exuy;
    uy(2:end-1,2:end-1) = Tuy(2:end-1,2:end-1);
    [ErrUyInf,ErrUy2] = cmpr(Exuy,uy,1);
    
    uxx = Exuxx;
    uxx(2:end-1,2:end-1) = Tuxx(2:end-1,2:end-1);
    [ErrUxxInf,ErrUxx2] = cmpr(Exuxx,uxx,1);
    
    uyy = Exuyy;
    uyy(2:end-1,2:end-1) = Tuyy(2:end-1,2:end-1);
    [ErrUyyInf,ErrUyy2] = cmpr(Exuyy,uyy,1);
    
    uxy = Exuxy;
    uxy(2:end-1,2:end-1) = Tuxy(2:end-1,2:end-1);
    [ErrUxyInf,ErrUxy2] = cmpr(Exuxy,uxy,1);
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
            
end
end


end
	

function [Linf,L2] = cmpr(ex,u,GG)

    if nargout==2, GG=0;end
    tmp = ex - u;
    
    if GG==1
        tmp(ExtPrb.Scatterer.GridGamma) = 0;
        u(ExtPrb.Scatterer.GridGamma)=0;
    end
    
    Linf = norm(tmp(:),inf)/norm(u(:),inf);
    L2   = norm(tmp(:),2)/norm(u(:),2);
end

function e = Exact(FocalDist,eta,phi,ExParams)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);
	
	e = (x.^ExParams.c) .* (y.^ExParams.d);
end

function dedx = dxdExact(FocalDist,eta,phi,ExParams)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);
	
	dedx = ExParams.c*(x.^(ExParams.c-1)) .* (y.^ExParams.d);
end

function dedy = dydExact(FocalDist,eta,phi,ExParams)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);
	
	dedy = ExParams.d*(x.^ExParams.c) .* (y.^(ExParams.d-1));
end

function dedn = dndExact(FocalDist,eta,phi,ExParams)
	%x  = FocalDist*cosh(eta).*cos(phi);
	%y  = FocalDist*sinh(eta).*sin(phi);
    xn = FocalDist*sinh(eta).*cos(phi);
    yn = FocalDist*cosh(eta).*sin(phi);
    
	dedn = dxdExact(FocalDist,eta,phi,ExParams).*xn +  dydExact(FocalDist,eta,phi,ExParams).*yn;
end
