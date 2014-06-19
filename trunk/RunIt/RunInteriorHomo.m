function RunInteriorHomo
    % Semyon Method, Homo HLM, Cart coord, 
    
%global R A x y x1 xn y1 yn dx dy cols rows NHR  ebinf etinf IntEta k0 k FocDist Eta%n 

%x1=-1.2;xn=1.2;
%y1=-1.2;yn=1.2;

 x1=-1.7;xn=1.2;
 y1=-1.7;yn=1.7;

Lx=xn-x1;Ly=yn-y1;
ebinf=[];etinf=[];

R0=1;
NHR=1.6;

a=1;%1.6;%2.5;
b=1/2;%0.8;

FocalDist = sqrt(a^2-b^2);
Eta0 = acosh(a/FocalDist);


%Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b));
%Parameterization  = Tools.Parameterizations.ParametricKite(struct('a',1,'b',.65*2,'c',1.5));
%Parameterization  = Tools.Parameterizations.ParametricSubmarine(struct('a',1,'b',1/2,'c',0,'p',200));
Parameterization  = Tools.Parameterizations.ParametricHeart(struct('a',13/16,'b',-5/16,'c',-2/16,'d',-1/16,'e',1,'p',3));

ScatType = 'StarShapedScatterer'; %'ellipse' or 'circle' or 'StarShapedScatterer'
BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
ChebyshevRange = struct('a',-pi,'b',pi);%don't change it

if strcmpi(ScatType,'ellipse')
    ExParams  = struct('ScattererType','ellipse','eta',Eta0,'FocalDistance',FocalDist);
elseif strcmpi(ScatType,'circle')
    ExParams  = struct('ScattererType','circle','r',R0);
elseif strcmpi(ScatType,'StarShapedScatterer')
    %Parameterization.XHandle = Tools.Parametirzations.AcosBtWD(a,1);
    %Parameterization.YHandle = Tools.Parametirzations.AsinBtWD(b,1);
    ExParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization);
end

fprintf('Method:%s,\t Ellipse: a=%d; \t b=%d \n',ScatType,a,b);

for k = 1%[1,5,10,15,20,25]
    
    ErrPre = 0; ErrXiPre =0;
    
    f   =@(th) Exact(th,k,ExParams);%(R0,th,k);
    dfdn=@(th) drExact(th,k,ExParams);%((R0,th,k);
	
	if strcmpi(BType,'Chebyshev')
		Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
	elseif strcmpi(BType,'Fourier')
		Basis = Tools.Basis.FourierBasis.BasisHelper(f,dfdn);
	end

for n=1:4 %run different grids
tic
	%build grid
% 		Nx=2.^(n+1)+5;	Ny=2.^(n+1)+5;
p=4;%3;
	Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
	%dx=Lx/(Nx-1); 	dy=Ly/(Ny-1);
    	
    Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
	
	% polar wavenumber isn't designed to work with ellipse scatterer
	% polar/elliptical wavenumber works with polar scatterer in general, but the problem here become inhomoginious, 
	% so generally to speak it is wrong place to use it here, mainly because of compariseon to exact here
    WaveNumberClsHandle = @Tools.Coeffs.ConstantWaveNumber;%WaveNumberElliptical;%ConstantWaveNumber;%WaveNumberPolarR
    WaveNumberAddParams = struct('k',k,'r0',NHR);
   
    if strcmpi(ScatType,'ellipse')
        ScattererClsHandle  = @Tools.Scatterer.EllipticScatterer;
        ScattererAddParams  = struct('Eta0',Eta0,'FocalDistance',FocalDist);
    elseif strcmpi(ScatType,'circle')
        ScattererClsHandle  = @Tools.Scatterer.PolarScatterer;
        ScattererAddParams  = struct('r0',R0,'ExpansionType',25);
    elseif strcmpi(ScatType,'StarShapedScatterer')
        ScattererClsHandle  = @Tools.Scatterer.StarShapedScatterer;
        ScattererAddParams  = ExParams;
    end

    CollectRhs = 1;
    
    IntPrb =  Solvers.InteriorHomoSolver ... 
        (Basis,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams,CollectRhs);
    
    Q0 = IntPrb.Q0;%(:,1:2*M+1);
    Q1 = IntPrb.Q1;%(:,2*M+2:4*M+2);
        
    cn1 =( Q1 \ ( -Q0*Basis.cn0)) ;
    
        
    xi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
    xi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1;
    ebinf(n)=0;

    if strcmpi(ScatType,'ellipse')
        ExParams2  = struct('ScattererType','ellipse','eta',IntPrb.Scatterer.eta,'FocalDistance',FocalDist);
        xiex = Exact(IntPrb.Scatterer.phi,k,ExParams2);%(FocalDist,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi,k0,NHR);
    elseif strcmpi(ScatType,'circle')
        ExParams2 =ExParams;
        ExParams2.r = IntPrb.Scatterer.r;
        xiex = Exact(IntPrb.Scatterer.th,k,ExParams2);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
    elseif strcmpi(ScatType,'StarShapedScatterer')
        
        ExParams2 = struct('ScattererType','StarShapedScatterer','r', IntPrb.Scatterer.r);        
        
        xiex = Exact(IntPrb.Scatterer.th,k,ExParams2);
        
    end
    
    %xiex = Exact(Eta(Split.GridGamma),phi,k0);
    ErrXi =norm(xiex -xi(IntPrb.GridGamma),inf);
    
    u = IntPrb.P_Omega(xi);
    
    t1=toc;
    
    % % % % % % % % % % % % % % % %
    % Comparison
    % % % % % % % % % % % % % % % %
    
    
    exact = zeros(Nx,Ny);   
    %exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.R(IntPrb.Scatterer.Np),IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k);
    if strcmpi(ScatType,'ellipse')
        ExParams3  = struct('ScattererType','ellipse','eta',IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),'FocalDistance',FocalDist);
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
        %TBD;
    elseif strcmpi(ScatType,'circle')
        ExParams3 =ExParams;
        ExParams3.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
     elseif strcmpi(ScatType,'StarShapedScatterer')
        
        ExParams3 = struct('ScattererType','StarShapedScatterer','r',IntPrb.Scatterer.R(IntPrb.Scatterer.Np));
        
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);        
    
    end
   
    
    
    %t2=toc;
    
    ErrTot =norm(exact(:)-u(:),inf);
    fprintf('k=%d,M=%d,N=%-4dx%-4d\t ErrXi=%d\t rate=%-4.2f ErrTot=%d\t rate=%-4.2f timeA=%d\n',...
        k,Basis.M, Nx,Ny,ErrXi,log2(ErrXiPre/ErrXi),ErrTot,log2(ErrPre/ErrTot),t1);
    ErrPre = ErrTot;
    ErrXiPre = ErrXi;
    
   % figure, plot(x,abs(u(fix(Ny/2),:)),'r+-',x,abs(exact(fix(Ny/2),:)),'bo');
%     saveas(gcf,['inhomo_k=' num2str(k) 'N=' num2str(Nx) '.jpg'],'jpg');
    % plot(x,abs(u(:,fix(Ny/2))),'r+-',x,abs(exact(:,fix(Ny/2))),'bo');

end

%Linf=log2(etinf(1:end-1)./etinf(2:end))
%Lbinf=log2(ebinf(1:end-1)./ebinf(2:end))
end
end



function g = Exact(th,k,Params)%(R,th,k)
if strcmpi(Params.ScattererType,'ellipse')
    x = Params.FocalDistance * cosh(Params.eta) .* cos(th);
    %          y = Params.FocalDistance * sinh(Params.eta) .* sin(th);
elseif strcmpi(Params.ScattererType,'circle')
    x = Params.r .* cos(th);
    %         y = Params.r .* sin(th);
elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
    try
        x = Params.Parameterization.XHandle.Derivatives(th);
    catch
        x = Params.r.*cos(th);
    end
end

%g = exp(1i.*k.*R.*cos(th));
g = exp(1i.*k.*x);
end


function dne = drExact(th,k,Params)%(R,th,k)
if strcmpi(Params.ScattererType,'ellipse')
    dx = Params.FocalDistance * sinh(Params.eta) .* cos(th);
    %         dy = Params.FocalDistance * cosh(Params.eta) .* sin(th);
elseif strcmpi(Params.ScattererType,'circle')
    dx = cos(th);
    %         dy = sin(th);
elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
    [x,dx] = Params.Parameterization.XHandle.Derivatives(th);
    [y,dy] = Params.Parameterization.YHandle.Derivatives(th);
end

e = Exact(th,k,Params);%(R,th,k);
%dne = 1i.*k.*cos(th).*e;
dne = 1i.*k.*dx.*e;


if strcmpi(Params.ScattererType,'StarShapedScatterer')
   h = sqrt(dx.^2 + dy.^2);
   % dne = dne./h;
	
	%fx dy + fy (-dx)
	
	dne  = (1i.*k.*e).*dy./h; %=fx dy, here fy=0
	
end


end
