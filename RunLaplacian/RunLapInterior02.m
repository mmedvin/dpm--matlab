function RunLapInterior02

    a=1;%2.5;
    b=1/2;	
	
	x1=-2;xn=2;  
    y1=-2;yn=2;
    
    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);

    BIn = 1e0;
    BOut = 1e3;

	BType		= 'Fourier';
    
     ExpansionType=35;
     Stencil = 13;
   
	f       =@(phi) IntExact(FocalDistance,Eta0,phi);
	dfdn    =@(phi) dndIntExact(FocalDistance,Eta0,phi);

	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,dfdn);
	end

for	   LinearSolverType = 1
    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end
    
    ErrUInfPre = 0; ErrU2Pre = 0; ErrUxInfPre = 0; ErrUx2Pre = 0; ErrUyInfPre = 0; ErrUy2Pre = 0; ErrUxxInfPre = 0; ErrUxx2Pre = 0; ErrUyyInfPre = 0; ErrUyy2Pre = 0; ErrUxyInfPre = 0; ErrUxy2Pre = 0;
	
    fprintf('Problem RunLapInterior02, M=%d, LinearSolverType = %d, BIn=%d, BOut=%d\n', Basis.M, LinearSolverType, BIn, BOut);
        
	for n=1:5 %run different grids
		tic
		%build grid
		p=4;%4;
		Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
        Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
        %------------------------------------------------------------------
		        		
		IntPrb = Solvers.InteriorHomoLaplacianSolver(struct(...
                 'Basis'             , Basis, ...
                 'Grid'              , Grid, ...
                 'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
                 'CoeffsParams'      , struct('a',BIn,'b',BIn,'sigma',0), ...
                 'ScattererHandle'   ,  @Tools.Scatterer.EllipticScatterer, ...
                 'ScattererParams'   , struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',ExpansionType, 'Stencil', Stencil), ...
                 'CollectRhs'        , CollectRhs, ...
                 'DiffOp'            , @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs, ...
                 'DiffOpParams'      , struct('BC_x1',0,'BC_xn', 0,'BC_y1',0,'BC_yn',0, 'LinearSolverType', LinearSolverType) ...
                 ));
		%------------------------------------------------------------------
	 
		
		%interior

		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		%xiex = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		
		% Neumann problem - not working need to add data for uniquiness
		%Intcn0 =( IntPrb.Q0 \ ( -IntPrb.Q1*Basis.cn1 )) ;
		%Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Intcn0 + IntPrb.W1(IntPrb.GridGamma,:)*Basis.cn1;
		%xiex(IntPrb.GridGamma) = IntExact(FocalDistance,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi);

		Intcn0 = Basis.cn0;
		IQ0 = IntPrb.Q0;
		IQ1 = IntPrb.Q1;
 		Intcn1 =( IQ1 \ ( -IQ0*Intcn0 )) ;    
 		Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Intcn0 + IntPrb.W1(IntPrb.GridGamma,:)*Intcn1;
    
		%xiex  = IntExact(FocalDistance,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi);
		%Intxi(ExtPrb.GridGamma) = xiex; %debig
		
		Intu = IntPrb.P_Omega(Intxi);

    %------------------------------------------------------------------
    
   t1=toc; 

   %------------------------------------------------------------------
    % Comparison
    %------------------------------------------------------------------
    
    Intexact = zeros(size(Grid.R));
	
    Intexact(IntPrb.Scatterer.Np) = IntExact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np));
    
    
    [ErrUInf,ErrU2] = cmpr(Intu(IntPrb.Scatterer.Mp),Intexact(IntPrb.Scatterer.Mp));
    
    CDInt = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy);

    [iTux,iTuy,iTuxx,iTuyy,iTuxy]       = CDInt.CartesianDerivatives(Intu);
    [iExux,iExuy,iExuxx,iExuyy,iExuxy]  = CDInt.CartesianDerivatives(Intexact);
    
    [ErrUxInf,ErrUx2]   = cmpr(iTux(IntPrb.Scatterer.Mp),iExux(IntPrb.Scatterer.Mp));
    [ErrUyInf,ErrUy2]   = cmpr(iTuy(IntPrb.Scatterer.Mp),iExuy(IntPrb.Scatterer.Mp));
	
    [ErrUxxInf,ErrUxx2] = cmpr(iTuxx(IntPrb.Scatterer.Mp),iExuxx(IntPrb.Scatterer.Mp));
    [ErrUyyInf,ErrUyy2] = cmpr(iTuyy(IntPrb.Scatterer.Mp),iExuyy(IntPrb.Scatterer.Mp));
    [ErrUxyInf,ErrUxy2] = cmpr(iTuxy(IntPrb.Scatterer.Mp),iExuxy(IntPrb.Scatterer.Mp));
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

    if nargin==2, GG=[];end
    tmp = ex - u;
    
    if ~isempty(GG)
        tmp(GG) = 0;
        u(GG)=0;
    end

    Linf = norm(tmp(:),inf);%/norm(u(:),inf);
    L2   = norm(tmp(:),2);%/norm(u(:),2);
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