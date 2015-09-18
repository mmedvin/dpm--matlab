function RunLapBL54
% 5.3  bealy-layton

    %a=0.9; b=0.7; %table 1
    a=0.9; b=0.1; %table 2
    
	x1=-1.1;xn=1.1;
	y1=-1.1;yn=1.1;

	Lx=xn-x1; Ly=yn-y1;

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
    
    Order=2;
    if Order==2 , ExpansionType=33; Stencil=5; else ExpansionType=35; 
        Stencil=9;%13;% 
    end
    
    ExParams.c=9;
    ExParams.d=8;

	BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
    ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
    
	f   =@(phi) IntExact(FocalDistance,Eta0,phi,ExParams);
	fn   =@(phi) dndIntExact(FocalDistance,Eta0,phi,ExParams);
	    
	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,fn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,fn,20);
	end

for	   LinearSolverType = 0
    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end

    ErrUInfPre = 0; ErrU2Pre = 0; ErrUxInfPre = 0; ErrUx2Pre = 0; ErrUyInfPre = 0; ErrUy2Pre = 0; 
    	
    fprintf('Problem BL54, M=%d, LinearSolverType = %d , Order=%d \n', Basis.M, LinearSolverType,Order);
        
	for n=1:4 %6 %run different grids
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

        Setup = struct( 'Basis'     , Basis, ...
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
                );
        
		IntPrb =  Solvers.InteriorLaplacianSolver(Setup);
							
        ExtPrb =  Solvers.ExteriorHomoLaplacianSolver(Setup);
		%------------------------------------------------------------------
	if 1
        nGGext = numel(ExtPrb.GridGamma);
        nGGint = numel(IntPrb.GridGamma);
        
		ZerosInt=spalloc(nGGint,Basis.NBss,0);
        ZerosExt=spalloc(nGGext,Basis.NBss,0);
        
		Eye = speye(Basis.NBss);
		Zeros2=spalloc(Basis.NBss,Basis.NBss,0);
		
		rhs = zeros(nGGext+nGGint + 2*Basis.NBss,1);
		rhs(1:nGGint)	= (-IntPrb.TrGF -IntPrb.Qf);
		%rhs((nGGint+1):(nGGext+nGGint))= -(ExtPrb.TrGF +ExtPrb.Qf);

		% [u_in] = x0^2 - y0^2 
		% [u_out] = sin(x0).* cos(y0) 
		% [ beta un_in] = 0
		% [beta un_out] = - 2 b sin(x0).* cos(y0)
		
		phi=linspace(0,2*pi, Basis.NBss+1);
		phi = phi(1:Basis.NBss);

        % 		x0 = FocalDistance*cosh(Eta0).*cos(phi);
        % 		y0 = FocalDistance*sinh(Eta0).*sin(phi);
        % 		xn0 = FocalDistance*sinh(Eta0).*cos(phi);
        % 		yn0 = FocalDistance*cosh(Eta0).*sin(phi);
		
        if 1
            ICu =  IntExact(FocalDistance,Eta0,phi,ExParams);
            ICuc = Tools.Basis.FourierBasis.FftCoefs(ICu, Basis.NBss);
            
            ICun =  dndIntExact(FocalDistance,Eta0,phi,ExParams);
            ICunc = Tools.Basis.FourierBasis.FftCoefs(ICun, Basis.NBss);
            
            rhs((nGGext+nGGint+1):(nGGext+nGGint+Basis.NBss))= ICuc;
            rhs((nGGext+nGGint+Basis.NBss+1):end)= ICunc;
            
            cn = [IntPrb.Q0,IntPrb.Q1,ZerosInt,ZerosInt;
                ZerosExt,ZerosExt,ExtPrb.Q0,ExtPrb.Q1;
                Eye,      Zeros2,   -Eye, Zeros2;
                Zeros2, Eye, Zeros2, -Eye]\rhs;
        else           
            ICuout = 0 ;
            ICuin= IntExact(FocalDistance,Eta0,phi,ExParams);
            ICuinc = Tools.Basis.FourierBasis.FftCoefs(ICuin, Basis.NBss);
            ICuoutc = Tools.Basis.FourierBasis.FftCoefs(ICuout, Basis.NBss);
            
            ICuoutn = 0 ;
            ICuinn = dndIntExact(FocalDistance,Eta0,phi,ExParams);
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
            
            %[Q,R,P]=qr(A'*A);
            %cn = P*( R\(Q\(A'*rhs)));
            cn=A\rhs;
        end
		Intcn=cn(1:2*Basis.NBss);
		Extcn= cn(2*Basis.NBss+1:end); 
		
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*Intcn + IntPrb.Wf(IntPrb.GridGamma);
		Intu = IntPrb.P_Omega(Intxi);
				
		Extxi = spalloc(Nx,Ny,length(ExtPrb.GridGamma));
		Extxi(ExtPrb.GridGamma) = ExtPrb.W(ExtPrb.GridGamma,:)*Extcn;%  + ExtPrb.Wf(ExtPrb.GridGamma));
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
		   
        t1=toc;

   %------------------------------------------------------------------
    % Comparison
    %------------------------------------------------------------------
    
   Intexact = zeros(size(Grid.R));
    Extexact = zeros(size(Grid.R));
	
    %Intu = zeros(size(Grid.R));
    
    Intexact(IntPrb.Scatterer.Np) = IntExact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np),ExParams);
    %Extexact(ExtPrb.Scatterer.Nm) = 0;% ExtExact(FocalDistance,ExtPrb.Scatterer.Eta(ExtPrb.Scatterer.Nm),ExtPrb.Scatterer.Phi(ExtPrb.Scatterer.Nm),ExParams);
    
    U=zeros(size(Grid.R));
    U(ExtPrb.Scatterer.Mm)=Extu(ExtPrb.Scatterer.Mm);
    U(IntPrb.Scatterer.Mp)=Intu(IntPrb.Scatterer.Mp);
    
    E=zeros(size(Grid.R));
    E(ExtPrb.Scatterer.Mm)=Extexact(ExtPrb.Scatterer.Mm);
    E(IntPrb.Scatterer.Mp)=Intexact(IntPrb.Scatterer.Mp);
    
    [ErrUInf,ErrU2] = cmpr(E,U,[]);
    
    CD = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy);
    [Ux,Uy]  = CD.CartesianDerivatives(U);    
    [Ex,Ey]  = CD.CartesianDerivatives(E);

    
    [ErrUxInf,ErrUx2]   = cmpr(Ex,Ux,[]);
    [ErrUyInf,ErrUy2]   = cmpr(Ey,Uy,[]);
	
%------------------------------------------------------------------
%      
%     fprintf('N=%-6dx%-7d Eu=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Eux=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euxx=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euyy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f Euxy=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f timeA=%-6.2f\n',...
%         Nx,Ny,ErrUInf,   ErrU2,   log2(ErrUInfPre/ErrUInf),     log2(ErrU2Pre/ErrU2), ...
%               ErrUxInf,  ErrUx2,  log2(ErrUxInfPre/ErrUxInf),   log2(ErrUx2Pre/ErrUx2), ...  
%               ErrUyInf,  ErrUy2,  log2(ErrUyInfPre/ErrUyInf),   log2(ErrUy2Pre/ErrUy2), ...  
%               ErrUxxInf, ErrUxx2, log2(ErrUxxInfPre/ErrUxxInf), log2(ErrUxx2Pre/ErrUxx2), ...
%               ErrUyyInf, ErrUyy2, log2(ErrUyyInfPre/ErrUyyInf), log2(ErrUyy2Pre/ErrUyy2), ...  
%               ErrUxyInf, ErrUxy2, log2(ErrUxyInfPre/ErrUxyInf), log2(ErrUxy2Pre/ErrUxy2), ...
%             t1);

%$16\times16$	&	$9.8821e-04$&	$-	 $ 	&	$3.3575e-04$&	$-	 $	&	$2.7225e-03$&	$-	 $ 	&	$3.6927e-05$&	$-	 $ 	&	$2.3633e-05$&	$-	 $	&	$4.6294e-05$&	$-	 $	\\

%  str = sprintf('$%-6d\\\\times%-7d$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ \\\\\\\\ \n',...
%        Nx-1,Ny-1,ErrUInf,  log2(ErrUInfPre/ErrUInf), ErrUxInf,log2(ErrUxInfPre/ErrUxInf), ErrUyInf,log2(ErrUyInfPre/ErrUyInf),...
%                    ErrU2,  log2(ErrU2Pre/ErrU2),     ErrUx2,   log2(ErrUx2Pre/ErrUx2), ErrUy2,  log2(ErrUy2Pre/ErrUy2) );

 str = sprintf('$%-6d\\\\times%-7d$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$ & $%-10.4d$ & $%-6.2f$  \\\\\\\\ \n',...
       Nx,Ny,ErrUInf,  log2(ErrUInfPre/ErrUInf), ErrUxInf,log2(ErrUxInfPre/ErrUxInf), ErrUyInf,log2(ErrUyInfPre/ErrUyInf) );

       ErrUInfPre = ErrUInf; ErrUxInfPre = ErrUxInf;  ErrUyInfPre = ErrUyInf; 
    
 fprintf(str);
 %fprintf(fileID,str);
            
end
end


end
	

function [Linf,L2] = cmpr(ex,u,GG)

    tmp = ex - u;
    
    if exist('GG','var')
        tmp(GG) = 0;
        u(GG)=0;    
    end
    
    Linf = norm(tmp(:),inf)/norm(u(:),inf);
    L2   = norm(tmp(:),2)/norm(u(:),2);
end

function e = IntExact(FocalDist,eta,phi,ExParams)
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

function dedn = dndIntExact(FocalDist,eta,phi,ExParams)
	%x  = FocalDist*cosh(eta).*cos(phi);
	%y  = FocalDist*sinh(eta).*sin(phi);
    xn = FocalDist*sinh(eta).*cos(phi);
    yn = FocalDist*cosh(eta).*sin(phi);
    
	dedn = dxdExact(FocalDist,eta,phi,ExParams).*xn +  dydExact(FocalDist,eta,phi,ExParams).*yn;
end
