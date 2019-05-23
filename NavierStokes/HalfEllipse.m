function HalfEllipse()
    x1=-1.2;xn=1.2;
    y1=0;yn=0.7;
    
    a=1;
    b=1/2;
    
    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
    
    ScatType = 'ellipse'; %'StarShapedScatterer';% 'StarShapedScatterer'; %'ellipse' or 'circle' or 'StarShapedScatterer'
    BType = 'Chebyshev';
    %ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
    ChebyshevRange = struct('a',0,'b',pi);%????
    
    fprintf('\n');
    fprintf('HalfEllipse, Grid: x1=%f, xn=%f, y1=%f, yn=%f \n ',x1,xn,y1,yn);
    LinearSolverType = 0;
    CollectRhs = 1;
    
     ErrPre = 0;
    ErrXiPre = 0;
      
    f   = @(phi) Exact(FocalDistance,Eta0,phi);
    fn  = @(phi) dndExact(FocalDistance,Eta0,phi);

    Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,fn,ChebyshevRange);
    
    for n=1:4 %run different grids
        tic
        p=4;%3;
        Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
        
        Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
        
        
        BC = log((Grid.x(1)-2).^2 + (Grid.y(2:end-1)-2).^2 );
         %BC = log((Grid.x(2:end-1)-2).^2 + (Grid.y(1)-2).^2 );
       % BC = log((Grid.x-2).^2 + (Grid.y(1)-2).^2 );
        
        Setup = struct( 'Basis'     , Basis, ...
            'Grid'              , Grid, ...
            'CoeffsHandle'      , @Tools.Coeffs.ConstLapCoeffs, ...
            'CoeffsParams'      , struct('a',1,'b',1,'sigma',0), ... %sigma is supposed to be negative here if not zero
            'ScattererHandle'   , @Tools.Scatterer.HalfEllipticScatterer, ...
            'ScattererParams'   , struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',33, 'Stencil', 5), ...
            'CollectRhs'        , CollectRhs, ...
            'DiffOp'            , @Tools.DifferentialOps.LaplacianOpBCinRhs, ...
            'DiffOpParams'      , struct(   'BC_y1',  BC, 'BC_yn',  0,'BC_x1',0 , 'BC_xn',0, 'LinearSolverType', LinearSolverType, 'Order',2), ...
            ... %'SourceHandle'      , @Tools.Source.LaplaceSource_BL54_Interior, ...
            ... %'SourceParams'      , ExParams, ...
            'Extension'         , @Tools.Extensions.TwoTupleExtension, ...
            'ExtensionParams',[] ...
            );
        
        %IntPrb =  Solvers.InteriorHomoLaplacianSolver(Setup);
        IntPrb =  Solvers.ExteriorHomoLaplacianSolver(Setup);
        
        Q0 = IntPrb.Q0;
        Q1 = IntPrb.Q1;
        %Qf = IntPrb.Qf;
        
        %TrGF = IntPrb.TrGF;
        %cn1 =( Q1 \ ( -Q0*Basis.cn0 - TrGF - Qf)) ;
        
        cn1 =( Q1 \ ( -Q0*Basis.cn0 )) ;
        
        xi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
        xi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1;% + IntPrb.Wf{1}(IntPrb.GridGamma);
        
        u = IntPrb.P_Omega(xi);
        
        %comparison
        
        exact = zeros(size(Grid.R));
        exact(IntPrb.Scatterer.Np) = Exact(FocalDistance,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np));
        
        ErrTot =norm(exact(:)-u(:),inf);
        fprintf('NBss0=%d,NBss1=%d,N=%-4dx%-4d\t ErrTot=%d\t rate=%-5.2f \n',...
            Basis.NBss0,Basis.NBss1, Nx,Ny,ErrTot,log2(ErrPre/ErrTot));
        ErrPre = ErrTot;
    end
end


function e = Exact(FocalDist,eta,phi)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);
	
	e = log((x-2).^2 + (y-2).^2 );	
end

function dedn = dndExact(FocalDist,eta,phi)
	x = FocalDist*cosh(eta).*cos(phi);
	y = FocalDist*sinh(eta).*sin(phi);

    dx = FocalDist * sinh(eta) .* cos(phi);
    dy = FocalDist * cosh(eta) .* sin(phi);
    
	dedn = (2*(x-2).*dx - 2*(y-2).*dy)./ ((x-2).^2 + (y-2).^2 );	
end

