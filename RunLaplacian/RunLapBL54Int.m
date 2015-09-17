function RunLapBL54
% 5.3  bealy-layton

    %a=0.9; b=0.7; %table 1
    a=0.9; b=0.1; %table 2
    
	x1=-1.1;xn=1.1;
	y1=-1.1;yn=1.1;

	Lx=xn-x1; Ly=yn-y1;

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
    
    Order=4;
    if Order==2 , ExpansionType=33; Stencil=5; else ExpansionType=35; 
        Stencil=9;%13;% 
    end
    
    ExParams.c=9;
    ExParams.d=8;

	BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
    ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
    
	f   =@(phi) Exact(FocalDistance,Eta0,phi,ExParams);
	fn   =@(phi) dndExact(FocalDistance,Eta0,phi,ExParams);
	    
	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,fn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,fn,20);
	end

for	   LinearSolverType = 0
    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end

    ErrUInfPre = 0; ErrU2Pre = 0; ErrUxInfPre = 0; ErrUx2Pre = 0; ErrUyInfPre = 0; ErrUy2Pre = 0; 
    	
    fprintf('Problem BL54, M=%d, LinearSolverType = %d , Order=%d \n', Basis.M, LinearSolverType,Order);
        
	for n=1:6 %6 %run different grids
		tic
		%build grid
        p=20;
		Nx=p*2^n; Ny=p*2^n;
        
        %p=5;
        %Nx=2.^(n+p);	Ny=2.^(n+p);
		
		Grid             = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);

        %------------------------------------------------------------------
        if Order==4
            %DiffOp = @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs;
            DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat4OrdrConst;
        elseif Order==2
            DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
            %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
        end

		%------------------------------------------------------------------
               		        
		IntPrb =  Solvers.InteriorLaplacianSolver(struct( ...
                  'Basis'             , Basis, ...
                  'Grid'              , Grid, ...
                  'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                  'CoeffsParams'      , struct('a',1,'b',1,'sigma',0), ...
                  'ScattererHandle'   , @Tools.Scatterer.EllipticScatterer, ...
                  'ScattererParams'   , struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',ExpansionType, 'Stencil', Stencil), ...
                  'CollectRhs'        , CollectRhs, ...
                  'DiffOp'            , DiffOp, ...
                  'DiffOpParams'      , struct(   'BC_y1',  0, 'BC_yn',  0,'BC_x1',0 , 'BC_xn',0, 'LinearSolverType', LinearSolverType, 'Order',Order), ...
                  'SourceHandle'      , @Tools.Source.LaplaceSource_BL54_Interior, ...
                  'SourceParams'      , ExParams ...
                  ));
		%------------------------------------------------------------------
		
		Q0 = IntPrb.Q0;
		Q1 = IntPrb.Q1;
		TrGF = IntPrb.TrGF;
		Qf = IntPrb.Qf;
		Intcn1 =( Q1 \ ( -Q0*Basis.cn0 -TrGF - Qf)) ;
		
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = ...
			IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*Intcn1 + IntPrb.Wf(IntPrb.GridGamma);
		
		%xiex  = Exact(FocalDistance,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi,ExParams);
		%Intxi(IntPrb.GridGamma) = xiex; %debug
		
		Extu = IntPrb.P_Omega(Intxi); 
		   
        t1=toc;

   %------------------------------------------------------------------
    % Comparison
    %------------------------------------------------------------------
    
    exact = zeros(size(Grid.R));
    exact(IntPrb.Scatterer.Np) = Exact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np),ExParams);
    %exact(IntPrb.Scatterer.Mp) = Exact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Mp),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Mp),ExParams);    
    
    u = zeros(size(Grid.R));
    u(IntPrb.Scatterer.Np) = Extu(IntPrb.Scatterer.Np);
    %u(IntPrb.Scatterer.Mp) = Extu(IntPrb.Scatterer.Mp);
    
	[ErrUInf,ErrU2] = cmpr(exact,u);   
 
    CD = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy);
	[ux,uy] = CD.CartesianDerivatives(u(:));
    [Exux,Exuy] = CD.CartesianDerivatives(exact(:));
    
    [ErrUxInf,ErrUx2] = cmpr(Exux(IntPrb.Scatterer.Np),ux(IntPrb.Scatterer.Np));   
    [ErrUyInf,ErrUy2] = cmpr(Exuy(IntPrb.Scatterer.Np),uy(IntPrb.Scatterer.Np));
    
    %[ErrUxInf,ErrUx2] = cmpr(Exux(IntPrb.Scatterer.Mp),ux(IntPrb.Scatterer.Mp));
    %[ErrUyInf,ErrUy2] = cmpr(Exuy(IntPrb.Scatterer.Mp),uy(IntPrb.Scatterer.Mp));
    
    %------------------------------------------------------------------
    
    fprintf('N=%-6dx%-7d Eu=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Eux=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f timeA=%-6.2f\n',...
        Nx,Ny,ErrUInf,   ErrU2,   log2(ErrUInfPre/ErrUInf),     log2(ErrU2Pre/ErrU2), ...
              ErrUxInf,  ErrUx2,  log2(ErrUxInfPre/ErrUxInf),   log2(ErrUx2Pre/ErrUx2), ...  
              ErrUyInf,  ErrUy2,  log2(ErrUyInfPre/ErrUyInf),   log2(ErrUy2Pre/ErrUy2), ...  
            t1);
    
    ErrUInfPre = ErrUInf; ErrU2Pre = ErrU2; ErrUxInfPre = ErrUxInf; ErrUx2Pre = ErrUx2; ErrUyInfPre = ErrUyInf; ErrUy2Pre = ErrUy2; 
            
end
end


end
	

function [Linf,L2] = cmpr(ex,u,GG)

    tmp = ex - u;
    
    if exist('GG','var')
        tmp(GG) = 0;
        u(GG)=0;    end
    
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
