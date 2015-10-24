function RunNested


x1=-2.2;xn=2.2;
y1=-2.2;yn=2.2;
%y1=-0.7;yn=0.7;

R0=1;
R1=2;

%kite
%x1=-1.7;xn=1.2;
%y1=-1.7;yn=1.7;

NHR = 1.6;%1.6;
% k=1;
a=1;
b=0.5;
%for b=0.1:0.1:0.9
    fprintf('Internal Inhomogeneous problem inside ellipse of AR=%d \n',a/b);
FocalDist = sqrt(a^2-b^2);
Eta0 =acosh(1/FocalDist);


%doesn't expected to work Parameterization  = Tools.Parameterizations.ParametricHeart(struct('a',13/16,'b',-5/16,'c',-2/16,'d',-1/16,'e',1,'p',3));
Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b));
%Parameterization  = Tools.Parameterizations.ParametricKite(struct('a',1,'b',.65*2,'c',1.5));
%Parameterization  = Tools.Parameterizations.ParametricSubmarine(struct('a',1,'b',1/5,'c',2,'p',100));
%Parameterization  = Tools.Parameterizations.ParametricStar();


ScatType = 'circle';%'StarShapedScatterer';%'StarShapedScatterer'; %'ellipse' or 'circle' or 'StarShapedScatterer'
BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
ChebyshevRange = struct('a',-pi,'b',pi);%don't change it

if strcmpi(ScatType,'ellipse')
    ExParams  = struct('ScattererType','ellipse','eta',Eta0,'FocalDistance',FocalDist,'r0',NHR);
elseif strcmpi(ScatType,'circle')
    ExParams{1}  = struct('ScattererType','circle', 'r', R0,'r0',NHR);
    ExParams{2}  = struct('ScattererType','circle', 'r', R1,'r0',NHR);
elseif strcmpi(ScatType,'StarShapedScatterer')
    ExParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization,'r0',NHR);
end

%fprintf('Method:%s,\t Ellipse: a=%d; \t b=%d \n',ScatType,a,b);
fprintf('Method:RunInterior-%s,\t using %s Basis, r0=%d,r1=%d \n',ScatType,BType,R0,R1);
fprintf('Grid: x1=%f, xn=%f, y1=%f, yn=%f \n %s \n',x1,xn,y1,yn, Parameterization.Print);

for k = 1%[1,5]% [1,3,5] %[1,5,10,15,20,25]

    ebinfPre=0; etinfPre=0; etinf = 0; ebinf = 0;
    
    f   =@(th,n) Exact  (th,k,ExParams{n});
    dfdn=@(th,n) drExact(th,k,ExParams{n});

	if strcmpi(BType,'Chebyshev')
		Basis = struct( 'Interior', Tools.Basis.ChebyshevBasis.BasisHelper(@(th) f(th,1),@(th) dfdn(th,1),ChebyshevRange), ...
                        'Exterior', Tools.Basis.ChebyshevBasis.BasisHelper(@(th) f(th,2),@(th) dfdn(th,2),ChebyshevRange) ...
                       );
    elseif strcmpi(BType,'Fourier')
		Basis = struct( 'Interior', Tools.Basis.FourierBasis.BasisHelper(@(th) f(th,1),@(th) dfdn(th,1)), ...
                        'Exterior', Tools.Basis.FourierBasis.BasisHelper(@(th) f(th,2),@(th) dfdn(th,2)) ...
                       );
	end
fprintf('NBss: %d,%d,%d,%d\n',Basis.Interior.NBss0,Basis.Interior.NBss1,Basis.Exterior.NBss0,Basis.Exterior.NBss1);

for n=1:4 %run different grids
tic
	%build grid
    p=4;
	Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;    
	    
    Grid             = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);

    if strcmpi(ScatType,'circle')
        WaveNumberHandle = @Tools.Coeffs.WaveNumberPolarR;
        ScattererHandle  = @Tools.Scatterer.NestedPolarScatterer;
        ScattererParams  = struct('r0',R0,'r1',R1, 'Stencil', 9);
        Source           = @Tools.Source.HelmholtzSourceR;
        SourceParams     =  [];
        Extension        = @Tools.Extensions.NestedExtension;
        ExtentionParams  = struct('IntExtension',@Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension,'ExtExtension',@Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension);
    elseif strcmpi(ScatType,'ellipse')
        WaveNumberHandle = @Tools.Coeffs.WaveNumberElliptical;
        ScattererHandle  = @Tools.Scatterer.EllipticScatterer;
        ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDist, 'Stencil', 9);
        Source           = @Tools.Source.HelmholtzSourceElps;%HelmholtzSourceElpsDescrete;%
        SourceParams     =  struct('Grid',Grid);  
        Extension = @Tools.Extensions.TwoTupleExtension;
    elseif strcmpi(ScatType,'StarShapedScatterer')
        WaveNumberHandle = @Tools.Coeffs.WaveNumberStarShaped;%WaveNumberStarShaped;%WaveNumberPolarR;
        Source           = @Tools.Source.HelmholtzSourceStarShaped;%HelmholtzSourceR;
        ScattererHandle  = @Tools.Scatterer.StarShapedScatterer;
        ScattererParams  = ExParams;
        ScattererParams.Stencil=9;
        SourceParams     = [];
        Extension        = @Tools.Extensions.TwoTupleExtension;
    end
   
        
    IntPrb = Solvers.InteriorSolver(struct( ...
            'Basis'             , Basis,...
            'Grid'              , Grid, ...
            'CoeffsHandle'      , WaveNumberHandle, ...
            'CoeffsParams'      , struct('k',k,'r0',NHR), ...
            'ScattererHandle'   , ScattererHandle, ...
            'ScattererParams'   , ScattererParams, ...
            'CollectRhs'        , 1, ... %i.e. yes
            'SourceHandle'      , Source, ...
            'SourceParams'      , SourceParams, ...
            'Extension'         , Extension, ...
            'ExtensionParams'   , ExtentionParams ...
            ));
        
     Cn	= [ IntPrb.Q{1,2}, IntPrb.Q{2,2}] \ ([ -IntPrb.Q{1,1},-IntPrb.Q{2,1}]*[Basis.Interior.cn0; Basis.Exterior.cn0]  - IntPrb.TrGF - IntPrb.Qf)  ;
     InteriorCn1 = Cn(1:Basis.Interior.NBss1);
     ExteriorCn1 = Cn((1+Basis.Interior.NBss1):end);

     xi = spalloc(Nx,Ny,length(IntPrb.GridGamma));     
     xi(IntPrb.GridGamma) = IntPrb.W{1,1}(IntPrb.GridGamma,:)*Basis.Interior.cn0 + IntPrb.W{1,2}(IntPrb.GridGamma,:)*InteriorCn1 ...
                          + IntPrb.W{2,1}(IntPrb.GridGamma,:)*Basis.Exterior.cn0 + IntPrb.W{2,2}(IntPrb.GridGamma,:)*ExteriorCn1 ...
                          + IntPrb.Wf(IntPrb.GridGamma);


     if strcmpi(ScatType,'circle')
         ExParams2 = ExParams;
         ExParams2{1}.r = IntPrb.Scatterer.r{1};
         ExParams2{2}.r = IntPrb.Scatterer.r{2};
         xiex1 = Exact(IntPrb.Scatterer.th{1},k,ExParams2{1});
         xiex2 = Exact(IntPrb.Scatterer.th{2},k,ExParams2{2});
         
         ErrXi1 =norm(xiex1 -xi(IntPrb.Scatterer.InteriorScatterer.GridGamma),inf);
         ErrXi2 =norm(xiex2 -xi(IntPrb.Scatterer.ExteriorScatterer.GridGamma),inf);
         %xi(IntPrb.Scatterrer.InteriorScatterer.GridGamma)=xiex1;
         %xi(IntPrb.Scatterrer.ExteriorScatterer.GridGamma)=xiex2;

    elseif strcmpi(ScatType,'ellipse')
        ExParams2 = ExParams;
        ExParams2.eta = IntPrb.Scatterer.eta;
        xiex = Exact(IntPrb.Scatterer.phi,k,ExParams2);
        ebinf =norm(xiex -xi(IntPrb.GridGamma),inf);
    elseif strcmpi(ScatType,'StarShapedScatterer')
        if 0
            ExParams2 = ExParams;
            ExParams2.r = IntPrb.Scatterer.r;
            xiex = Exact(IntPrb.Scatterer.th,k,ExParams2);
            ebinf =norm(xiex -xi(IntPrb.GridGamma),inf);
        else
            if n > 1
                xiPre2= spalloc(Nx,Ny,nnz(xiPre));
                xiPre2(IntPrb.GridGamma) = xiPre(IntPrb.GridGamma);
                
                tmp = xi(1:2:end,1:2:end)-xiPre2(1:2:end,1:2:end);
                ebinf =norm(tmp(:),inf);
            end
            xiPre=spalloc(Nx*2-1,Ny*2-1,nnz(xi));
            xiPre(1:2:end,1:2:end) = xi;
        end
    end
    
    %xi(IntPrb.GridGamma)=xiex;
    
    
    u = IntPrb.P_Omega(xi);
    
    t1=toc;
    
    % % % % % % % % % % % % % % % %
    % Comparison
    % % % % % % % % % % % % % % % %
    
    
    exact = zeros(Nx,Ny);
     if strcmpi(ScatType,'circle')
        ExParams3 = ExParams;
        ExParams3.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
        etinf = norm(exact(:)-u(:),inf);
    elseif strcmpi(ScatType,'ellipse')
        ExParams3 = ExParams;
        ExParams3.eta = IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np);
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
        etinf = norm(exact(:)-u(:),inf);
    elseif strcmpi(ScatType,'StarShapedScatterer')
        if 0
            ExParams3 = ExParams;
            ExParams3.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
            exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);
            etinf = norm(exact(:)-u(:),inf);
        else
            if n > 1
                u1= spalloc(Nx,Ny,nnz(u0));
                u1(IntPrb.Np) = u0(IntPrb.Np);
                
                tmp = u(1:2:end,1:2:end)-u1(1:2:end,1:2:end);
                etinf =norm(tmp(:),inf);
            end
            
            u0=spalloc(Nx*2-1,Ny*2-1,nnz(u));
            u0(1:2:end,1:2:end)=u;
        end
    end
    t2=toc;
    
    
        fprintf('k=%d,NBss0=%d,NBss1=%d,N=%-10dx%-10d\t ErrXi=%d|%d\t rate=%-5.2f|%-5.2f etinf=%d|%-5.4f\ttimeA=%d\ttimeE=%d\n', ...
                 k,Basis.NBss0,Basis.NBss1, Nx,Ny,ErrXi1,ErrXi2,log2(ErrXi1Pre/ErrXi1),log2(ErrXi2Pre/ErrXi2), etinf,log2(etinfPre./etinf),t1,t2-t1);
    etinfPre = etinf;
    ErrXi1Pre = ErrXi1;
    ErrXi2Pre = ErrXi2;
    
end

fprintf('\n')
end
%end
end

function g = Exact(th,k0,Params)
if strcmpi(Params.ScattererType,'ellipse')
    k = Tools.Coeffs.WaveNumberElliptical.kkn(Params.FocalDistance,Params.eta,th,k0,Params.r0);
    x = Params.FocalDistance * cosh(Params.eta) .* cos(th);
    %          y = Params.FocalDistance * sinh(Params.eta) .* sin(th);
elseif strcmpi(Params.ScattererType,'circle')
    x = Params.r .* cos(th);
    %         y = Params.r .* sin(th);
     k = Tools.Coeffs.WaveNumberPolarR.kkr(Params.r,Params.r0,k0);
elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
    try
        x = Params.Parameterization.XHandle.Derivatives(th);
        y = Params.Parameterization.YHandle.Derivatives(th);
    catch
        x = Params.r.*cos(th);
        y = Params.r.*sin(th);
    end
    r=sqrt(x.^2+y.^2);
    k = Tools.Coeffs.WaveNumberPolarR.kkr(r,Params.r0,k0);
    
end

%g = exp(1i.*k.*R.*cos(th));
g = exp(1i.*k.*x);
end


function dne = drExact(th,k0,Params)%(R,th,k)
if strcmpi(Params.ScattererType,'ellipse')
    [k,kn] = Tools.Coeffs.WaveNumberElliptical.kkn(Params.FocalDistance,Params.eta,th,k0,Params.r0);
    kx=kn;
    x = Params.FocalDistance * cosh(Params.eta) .* cos(th);
    dx = Params.FocalDistance * sinh(Params.eta) .* cos(th);
    %         dy = Params.FocalDistance * cosh(Params.eta) .* sin(th);
elseif strcmpi(Params.ScattererType,'circle')
    x = Params.r .* cos(th);
    dx = cos(th);
    %         dy = sin(th);
    [k,kr] = Tools.Coeffs.WaveNumberPolarR.kkr(Params.r,Params.r0,k0);
    kx=kr;
elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
    [x,dx] = Params.Parameterization.XHandle.Derivatives(th);
    [y,dy] = Params.Parameterization.YHandle.Derivatives(th);
    
    r=sqrt(x.^2+y.^2);
    [k,kr] = Tools.Coeffs.WaveNumberPolarR.kkr(r,Params.r0,k0);        
    kx = kr.*x./r;
    %ky = kr.*y./r;
    
end

e = Exact(th,k,Params);%(R,th,k);
%dne = 1i.*k.*cos(th).*e;
dne = 1i.*(k.*dx+kx.*x).*e;


%dne = 1i.*FocDist.*cos(phi).*(kn.*cosh(eta)+k.*sinh(eta)).*e;

if strcmpi(Params.ScattererType,'StarShapedScatterer')
   h = sqrt(dx.^2 + dy.^2);
   % dne = dne./h;
	
	%fx dy + fy (-dx)
	
	dne  = 1i.*e.*(k+x.^2.*kr./r).*dy./h - (1i.*e.*x.*y.*kr./r).*dx./h;
	
end


end