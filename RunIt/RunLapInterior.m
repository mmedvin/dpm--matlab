function RunLapInterior
    % DPM for Non Homo (Lap-sigma)u=f, Cart coord
    
x1=-1;xn=1;
y1=-1;yn=1;
Lx=xn-x1;Ly=yn-y1;
ebinf=[];etinf=[];



ExParams.B  = 10^-3;
ExParams.C  = 0.1;
ExParams.r0 = 1/4;
ExParams.r  = ExParams.r0;

    
    f   =@(theta) Exact(ExParams,theta);
    dfdn=@(theta) drExact(ExParams,theta);

    Basis = Tools.Basis.FourierBasis.BasisHelper(f,dfdn);

for n=1:4 %run different grids
tic
	
    p=3;%3;
	Nx=2.^(n+p);	Ny=2.^(n+p);
    
    Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);

    CoeffsClsHandle = @Tools.Coeffs.LaplaceCoeffsPolar;
    CoeffsClsAddParams = ExParams;
    
    ScattererClsHandle  = @Tools.Scatterer.PolarScatterer;
    ScattererAddParams  = struct('r0',ExParams.r0,'ExpansionType',33);
    
    Source              = @Tools.Source.LaplaceSource;
    
    IntPrb = Solvers.InteriorSolver ... 
        (Basis,Grid,CoeffsClsHandle,CoeffsClsAddParams,ScattererClsHandle,ScattererAddParams,Source);
    
    Q0 = IntPrb.Q0;
    Q1 = IntPrb.Q1;
    Qf = IntPrb.Qf;
    
    TrGF = IntPrb.TrGF;
    
    cn1 =( Q1 \ ( -Q0*Basis.cn0 - TrGF - Qf)) ;
    
        
    xi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
    xi(IntPrb.GridGamma) = ...
        IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1 + IntPrb.Wf(IntPrb.GridGamma);
    
    ebinf(n)=0;
    xiex = Exact(FocalDist,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi,k0,NHR);
    %xiex = Exact(Eta(Split.GridGamma),phi,k0);
    ebinf(n) =norm(xiex -xi(IntPrb.GridGamma),inf);

    u = IntPrb.P_Omega(xi);
    
    %%%%%%%%%%%%%%
    
   t1=toc; 

   % % % % % % % % % % % % % % % %
    % Comparison
    % % % % % % % % % % % % % % % %
    
    exact = zeros(Nx,Ny);
        
    TestParams = ExParams;
    TestParams = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
    
    exact(IntPrb.Scatterer.Np) = Exact(TestParams,IntPrb.Scatterer.Th(IntPrb.Scatterer.Np));
       
    t2=toc;
    
    etinf(n) =norm(exact(:)-u(:),inf);
    fprintf('coeffs=%d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',0,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
        
end

Linf=log2(etinf(1:end-1)./etinf(2:end))
Lbinf=log2(ebinf(1:end-1)./ebinf(2:end))

%end
end


function e = Exact(Params,theta)
    r0 = 2*Params.r0;
    r  = Params.r   ;
    
    e   = ones(size(theta)).*r.^2;
    if length(r)==1 && r>r0
        e   = (1 - 1/8/Params.B - 1/Params.B)/4 + ( (r.^4)/2 + r.^2 )/Params.B + Params.C*log(2*r)/Params.B;        
    else
        e(r>r0)     = (1 - 1/8/Params.B - 1/Params.B)/4 + ( (r(r>r0).^4)/2 + r(r>r0).^2 )/Params.B + Params.C*log(2*r(r>r0))/Params.B;
    end
end

function er = drExact(Params,theta)
    r0 = 2*Params.r0;
    r  = Params.r   ;
  
    er  = ones(size(theta)).*2*r;
    
    if length(r)==1 && r>r0    
        er  = 2*( r.^3 + r )/Params.B + Params.C/Params.B./r;
    else
        er(r>r0)    = 2*( r(r>r0).^3 + r(r>r0) )/Params.B + Params.C/Params.B./r(r>r0);
    end
end
