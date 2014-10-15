function RunLaplacian01

    a=0.9;%2.5;
    b=0.1;	
	
    %a=0.9; b=0.1;
	
	x1=-1.1;xn=1.1;
	y1i=-0.5;yni=0.5;
	y1e=-1.1;yne=1.1;
	%Lx=xn-x1;Ly=yn-y1;    

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);

%        Coeff.a0 = 10;
%        Coeff.b0 = 10;
%        Coeff.a1 = 1;
%        Coeff.b1 = 1;
%     Coeff.c0 = 2;
%     Coeff.d0 = 1;
%      Coeff.c1 = 1;
    %Coeff.d1 = 1;

    ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
    
	BType		= 'Fourier'; %'Chebyshev';%'Fourier';
    
    %ExParams.CoeffsOut	= Tools.Coeffs.ConstLapCoeffs([], struct('a',Coeff.c0,'b',Coeff.d0,'sigma',0) ); 
    ExParams.CoeffsOut	= Tools.Coeffs.ExactCoeffs( struct('a',1,'b',2,'c',0,'d',3));
    ExParams.CoeffsIn	= Tools.Coeffs.ExactCoeffs( struct('a',2,'b',1,'c',0,'d',2));

%ExParams.CoeffsOut	= Tools.Coeffs.ConstLapCoeffs([], struct('a',Coeff.c0,'b',Coeff.d0,'sigma',0) );    
%ExParams.CoeffsIn	= Tools.Coeffs.ConstLapCoeffs([], struct('a',Coeff.c1,'b',Coeff.d1,'sigma',0) );

    ExParams.Eta0       = Eta0;
    
    f   =@(phi) ExtExact(FocalDistance,Eta0,phi,ExParams);
    fn  =@(phi) dndExtExact(FocalDistance,Eta0,phi,ExParams);    
	g   =@(phi) IntExact(FocalDistance,Eta0,phi,ExParams);
	gn  =@(phi) dndIntExact(FocalDistance,Eta0,phi,ExParams);
	
    

	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,g,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
        
        Basis0 = Tools.Basis.FourierBasis.BasisHelper(fn,gn);
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,g);
        if Basis0.NBss > Basis.NBss
            Basis = Tools.Basis.FourierBasis.BasisHelper(f,g,Basis0.NBss);
        end
        clear Basis0;
	end

for	   LinearSolverType = 0
    if LinearSolverType==0, CollectRhs = 1; else CollectRhs = 0;  end

    %ErrIntPre = 0; 	ErrExtPre = 0;	ErrTotPre = 0; 
    ErrUInfPre = 0; ErrU2Pre = 0; ErrUxInfPre = 0; ErrUx2Pre = 0; ErrUyInfPre = 0; ErrUy2Pre = 0; ErrUxxInfPre = 0; ErrUxx2Pre = 0; ErrUyyInfPre = 0; ErrUyy2Pre = 0; ErrUxyInfPre = 0; ErrUxy2Pre = 0;

    fprintf('Problem01, M=%d, LinearSolverType = %d, elps_a=%d, elps_b=%d\n', Basis.M, LinearSolverType, a, b);
    
    
	for n=1:4 %run different grids
		tic
		%build grid
		%p=3;%3;
		%Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
        %Nx=2.^(n+p);	Ny=2.^(n+p);
		
        p=20;
		Nx=p*2^n; Ny=p*2^n;        
        
		GridExt                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1e,yne,Ny);
		GridInt                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1i,yni,Ny);
		ScattererHandle  = @Tools.Scatterer.EllipticScatterer;
		ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',33);

		%------------------------------------------------------------------
		%InteriorCoeffsHandle = @Tools.Coeffs.ConstLapCoeffs;
		%InteriorCoeffsParams = struct('a',Coeff.a1,'b',Coeff.b1,'sigma',0);
        
        InteriorCoeffsHandle = @Tools.Coeffs.LaplaceCoeffsEllps1;
        InteriorCoeffsParams = struct('ca',2,'da',1,'ea',3,'WithB',1,'cb',4,'db',3,'eb',2,'sigma',0);
        
        InteriorSource          = @Tools.Source.LaplaceSource01_Interior;
        
        %------------------------------------------------------------------

		SourceParams = ExParams;
        		
        DiffOp = @Tools.DifferentialOps.LaplacianOpBCinRhs;
        %DiffOp = @Tools.DifferentialOps.LaplacianOpBCinMat;
                        
        DiffOpParamsInt = struct('BC_x1',0,'BC_xn', 0,'BC_y1',0,'BC_yn',0, 'LinearSolverType', LinearSolverType);
        
		IntPrb = Solvers.InteriorLaplacianSolver ...
			(Basis,GridInt,InteriorCoeffsHandle,InteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,InteriorSource,SourceParams,DiffOp,DiffOpParamsInt);
        
        %------------------------------------------------------------------
        %eZ = acosh((GridExt.X + 1i*GridExt.Y)/FocalDistance);
        %BeZ = [eZ(2:end-1,1), eZ(2:end-1,GridExt.Ny), eZ(1,2:end-1).' , eZ(GridExt.Nx,2:end-1).' ];
        %------BeZ = [eZ(1,2:end-1), eZ(GridExt.Ny,2:end-1), eZ(2:end-1,1).' , eZ(2:end-1,GridExt.Nx).' ];
        [BEta,BPhi] = GridExt.ToElliptical(FocalDistance);

        Boundaries.FocalDistance = FocalDistance;
        Boundaries.Eta0          = Eta0;
        Boundaries.Eta           = [BEta(2:end-1,1), BEta(2:end-1,end), BEta(1,2:end-1).' , BEta(end,2:end-1).' ];%real(BeZ);
        Boundaries.Phi           = [BPhi(2:end-1,1), BPhi(2:end-1,end), BPhi(1,2:end-1).' , BPhi(end,2:end-1).' ];%imag(BeZ);

        %Boundaries.Eta           = [BEta(1,2:end-1), BEta(end,2:end-1), BEta(2:end-1,1).' , BEta(2:end-1,end).' ];
        %Boundaries.Phi           = [BPhi(1,2:end-1), BPhi(end,2:end-1), BPhi(2:end-1,1).' , BPhi(2:end-1,end).' ];
        
%        Boundaries.Eta           = real(BeZ);
%        Boundaries.Phi           = imag(BeZ);
        
        %Exact.u	= ExtExact(FocalDistance,Boundaries.Eta,Boundaries.Phi,ExParams);
        Exact	= Tools.Exact.ExLapElps01(Boundaries, ExParams);
        DiffOpParamsExt = struct('BC_y1', Exact.u(:,1),'BC_yn', Exact.u(:,2),'BC_x1',Exact.u(:,3).','BC_xn',Exact.u(:,4).', 'LinearSolverType', LinearSolverType);
        
		        
%         c = ExParams.CoeffsOut.Derivatives('c',GridExt.x,GridExt.y);
%         d = ExParams.CoeffsOut.Derivatives('d',GridExt.x,GridExt.y);
%         
%         DiffOpParamsExt2 = struct('BC_y1', sin(c.*GridExt.x(1))      .* sin(d.*GridExt.y(2:end-1)).', ...
%                                  'BC_yn', sin(c.*GridExt.x(end))    .* sin(d.*GridExt.y(2:end-1)).',...
%                                  'BC_x1',(sin(c.*GridExt.x(2:end-1)).* sin(d.*GridExt.y(1))), ...
%                                  'BC_xn',(sin(c.*GridExt.x(2:end-1)).* sin(d.*GridExt.y(end))), 'LinearSolverType', LinearSolverType);
%         
        
        
        %ExteriorCoeffsHandle = @Tools.Coeffs.ConstLapCoeffs;
        %ExteriorCoeffsParams = struct('a',Coeff.a0,'b',Coeff.b0,'sigma',0);
        
        ExteriorCoeffsHandle = @Tools.Coeffs.LaplaceCoeffsEllps1;
        ExteriorCoeffsParams = struct('ca',5,'da',1,'ea',1000000,'WithB',1,'cb',1,'db',3,'eb',1000000,'sigma',0);       
        
        ExteriorSource       = @Tools.Source.LaplaceSource01_Exterior;
       
        ExtPrb =  Solvers.ExteriorLaplacianSolver ...
			(Basis,GridExt,ExteriorCoeffsHandle,ExteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,ExteriorSource,SourceParams,DiffOp,DiffOpParamsExt);
							
		%------------------------------------------------------------------
	if 1
		Zeros1=spalloc(numel(IntPrb.GridGamma),Basis.NBss,0);
        Zeros2=spalloc(numel(ExtPrb.GridGamma),Basis.NBss,0);
		Eye = speye(Basis.NBss);
		Zeros3=spalloc(Basis.NBss,Basis.NBss,0);
		
		nGGInt = numel(IntPrb.GridGamma);
        nGGExt = numel(ExtPrb.GridGamma);
		rhs = zeros(nGGInt + nGGExt + 2*Basis.NBss,1);
		rhs(1:nGGInt)	= (-IntPrb.TrGF -IntPrb.Qf);
		rhs(nGGInt+1:(nGGInt+nGGExt))= (-ExtPrb.TrGF -ExtPrb.Qf);

        Exact	= Tools.Exact.ExLapElps01(IntPrb.Scatterer, ExParams);
		
		[ICu,ICun] = Exact.InterfaceCondition(Basis.NBss); 
		ICuc = Tools.Basis.FourierBasis.FftCoefs(ICu, Basis.NBss);
		ICunc = Tools.Basis.FourierBasis.FftCoefs(ICun, Basis.NBss);
		
		rhs((nGGInt + nGGExt+1):(nGGInt + nGGExt+Basis.NBss))= ICuc;
		rhs((nGGInt + nGGExt+Basis.NBss+1):end)= ICunc;
		
 		cn = [IntPrb.Q0, IntPrb.Q1, Zeros1,      Zeros1    ;   
			  Zeros2,     Zeros2,     ExtPrb.Q0,  ExtPrb.Q1; 
			  Eye,       Zeros3,   -Eye,        Zeros3; 
			  Zeros3,    Eye,       Zeros3,    -Eye]\rhs;

		Intcn=cn(1:2*Basis.NBss);
		Extcn= cn(2*Basis.NBss+1:end); 
		
		Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
		Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*Intcn + IntPrb.Wf(IntPrb.GridGamma);
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
		
		%xiex  = ExtExact(FocalDistance,ExtPrb.Scatterer.eta,ExtPrb.Scatterer.phi,ExParams);
		%Extxi(ExtPrb.GridGamma) = xiex; %debug
		
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
		ITrGF = IntPrb.TrGF;
		IQf   = IntPrb.Qf;
 		Intcn1 =( IQ1 \ ( -IQ0*Intcn0 -ITrGF - IQf)) ;    
 		Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Intcn0 + IntPrb.W1(IntPrb.GridGamma,:)*Intcn1  + IntPrb.Wf(IntPrb.GridGamma);
    
		%xiex  = IntExact(FocalDistance,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi,ExParams);
		%Intxi(IntPrb.GridGamma) = xiex; %debug
		
		Intu = IntPrb.P_Omega(Intxi);
	end

    
   t1=toc; 

   %------------------------------------------------------------------
   % Comparison
   %------------------------------------------------------------------
   
   Intexact = zeros(size(GridInt.R));
   Extexact = zeros(size(GridExt.R));
   
   
   Intexact(IntPrb.Scatterer.Np) = IntExact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np),ExParams);
   Extexact(ExtPrb.Scatterer.Nm) = ExtExact(FocalDistance,ExtPrb.Scatterer.Eta(ExtPrb.Scatterer.Nm),ExtPrb.Scatterer.Phi(ExtPrb.Scatterer.Nm),ExParams);

    Extu(1,:)   = Extexact(1,:);
    Extu(end,:) = Extexact(end,:);
    Extu(:,1)   = Extexact(:,1);
    Extu(:,end) = Extexact(:,end);
   
    if 1
    [ErrUInf,ErrU2] = cmpr([Intu(IntPrb.Scatterer.Mp);Extu(ExtPrb.Scatterer.Mm)],[Intexact(IntPrb.Scatterer.Mp);Extexact(ExtPrb.Scatterer.Mm)],0);

    CDExt = Tools.Common.SecondDerivative(GridExt.Nx,GridExt.Ny,GridExt.dx,GridExt.dy);
    CDInt = Tools.Common.SecondDerivative(GridInt.Nx,GridInt.Ny,GridInt.dx,GridInt.dy);
    
    [eTux,eTuy,eTuxx,eTuyy,eTuxy]       = CDExt.CartesianDerivatives(Extu);
    [iTux,iTuy,iTuxx,iTuyy,iTuxy]       = CDInt.CartesianDerivatives(Intu);
    [eExux,eExuy,eExuxx,eExuyy,eExuxy]  = CDExt.CartesianDerivatives(Extexact);
    [iExux,iExuy,iExuxx,iExuyy,iExuxy]  = CDInt.CartesianDerivatives(Intexact);
            
    [ErrUxInf,ErrUx2]   = cmpr([iTux(IntPrb.Scatterer.Mp);eTux(ExtPrb.Scatterer.Mm)],[iExux(IntPrb.Scatterer.Mp);eExux(ExtPrb.Scatterer.Mm)],0);
    [ErrUyInf,ErrUy2]   = cmpr([iTuy(IntPrb.Scatterer.Mp);eTuy(ExtPrb.Scatterer.Mm)],[iExuy(IntPrb.Scatterer.Mp);eExuy(ExtPrb.Scatterer.Mm)],0);
    
    [ErrUxxInf,ErrUxx2] = cmpr([iTuxx(IntPrb.Scatterer.Mp);eTuxx(ExtPrb.Scatterer.Mm)],[iExuxx(IntPrb.Scatterer.Mp);eExuxx(ExtPrb.Scatterer.Mm)],0);
    [ErrUyyInf,ErrUyy2] = cmpr([iTuyy(IntPrb.Scatterer.Mp);eTuyy(ExtPrb.Scatterer.Mm)],[iExuyy(IntPrb.Scatterer.Mp);eExuyy(ExtPrb.Scatterer.Mm)],0);
    [ErrUxyInf,ErrUxy2] = cmpr([iTuxy(IntPrb.Scatterer.Mp);eTuxy(ExtPrb.Scatterer.Mm)],[iExuxy(IntPrb.Scatterer.Mp);eExuxy(ExtPrb.Scatterer.Mm)],0);
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
    else

         [ErrUInf,ErrU2] = cmpr(Intu(IntPrb.Scatterer.Mp),Intexact(IntPrb.Scatterer.Mp),0);
         [ErrUxInf,ErrUx2] = cmpr(Extu(ExtPrb.Scatterer.Mm),Extexact(ExtPrb.Scatterer.Mm),1);
         %[ErrUxInf,ErrUx2] = cmpr(Extu(2:end-1,2:end-1),Extexact(2:end-1,2:end-1),1);

          fprintf('N=%-6dx%-7d ErrInt=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f ErrExt=%-10.4d|%-10.4d rt=%-6.2f|%-6.2f timeA=%-6.2f\n',...
            Nx,Ny,ErrUInf,   ErrU2,   log2(ErrUInfPre/ErrUInf),     log2(ErrU2Pre/ErrU2), ...
            ErrUxInf,  ErrUx2,  log2(ErrUxInfPre/ErrUxInf),   log2(ErrUx2Pre/ErrUx2), ...
            t1);
        
        ErrUInfPre = ErrUInf; ErrU2Pre = ErrU2; ErrUxInfPre = ErrUxInf; ErrUx2Pre = ErrUx2;
    end
    
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


function e = IntExact(FocalDist,eta,phi,ExParams)

	x  = FocalDist*cosh(eta).*cos(phi);
	y  = FocalDist*sinh(eta).*sin(phi);
	
    c = ExParams.CoeffsIn.Derivatives('c',x,y);
    d = ExParams.CoeffsIn.Derivatives('d',x,y);

    e = sin(c.*x).*sin(d.*y);
    
end
function dnde = dndIntExact(FocalDist,eta,phi,ExParams)
	x  = FocalDist*cosh(eta).*cos(phi);
	y  = FocalDist*sinh(eta).*sin(phi);
    xn = FocalDist*sinh(eta).*cos(phi);
    yn = FocalDist*cosh(eta).*sin(phi);

    [c,cy] = ExParams.CoeffsIn.Derivatives('c',x,y);
    [d,dx] = ExParams.CoeffsIn.Derivatives('d',x,y);
    

    ux = y.*dx.*sin(c.*x).*cos(d.*y) + c.*cos(c.*x).*sin(d.*y);
    uy = x.*cy.*cos(c.*x).*sin(d.*y) + d.*cos(d.*y).*sin(c.*x);
    
    dnde = ux.*xn + uy.*yn;  
end

function dnde = dndExtExact(FocalDist,eta,phi,ExParams)
	x  = FocalDist*cosh(eta).*cos(phi);
	y  = FocalDist*sinh(eta).*sin(phi);
    xn = FocalDist*sinh(eta).*cos(phi);
    yn = FocalDist*cosh(eta).*sin(phi);

    [c,cy] = ExParams.CoeffsOut.Derivatives('c',x,y);
    [d,dx] = ExParams.CoeffsOut.Derivatives('d',x,y);
    
    ux = y.*dx.*sin(c.*x).*cos(d.*y) + c.*cos(c.*x).*sin(d.*y);
    uy = x.*cy.*cos(c.*x).*sin(d.*y) + d.*cos(d.*y).*sin(c.*x);
    
    dnde = ux.*xn + uy.*yn;  
end

function e = ExtExact(FocalDist,eta,phi,ExParams)
	x  = FocalDist*cosh(eta).*cos(phi);
	y  = FocalDist*sinh(eta).*sin(phi);

    c = ExParams.CoeffsOut.Derivatives('c',x,y);
    d = ExParams.CoeffsOut.Derivatives('d',x,y);
    
    e = sin(c.*x).*sin(d.*y);
end


