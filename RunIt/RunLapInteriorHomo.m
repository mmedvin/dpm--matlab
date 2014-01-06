function RunLapInteriorHomo
    % Semyon Method, Homo HLM, Cart coord, 
    
%global R A x y x1 xn y1 yn dx dy cols rows NHR  ebinf etinf IntEta k0 k FocDist Eta%n 

x1=-1.2;xn=1.2;
y1=-1.2;yn=1.2;
Lx=xn-x1;Ly=yn-y1;
ebinf=[];etinf=[];

R0=1;

a=1;%1.6;%2.5;
b=1/2;%0.8;
ellipse = struct('a',a,'b',b);
tower = struct('c',0,'p',20);
FocalDist = sqrt(a^2-b^2);
Eta0 = acosh(a/FocalDist);


ScatType = 'circle'; %'ellipse' or 'circle' or 'submarine'
BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
ChebyshevRange = struct('a',-pi,'b',pi);%don't change it

LapCoeffs = struct('a',1,'b',1,'sigma',-2);

if strcmpi(ScatType,'ellipse')
    ExParams  = struct('ScattererType','ellipse','eta',Eta0,'FocalDistance',FocalDist);
elseif strcmpi(ScatType,'circle')
    ExParams  = struct('ScattererType','circle','r',R0);
elseif strcmpi(ScatType,'submarine')
    ExParams  = struct('ScattererType','submarine','ellipse',ellipse,'tower',tower);
end

    
    f   =@(th) Exact(th,LapCoeffs,ExParams);%(R0,th,k);
    dfdn=@(th) drExact(th,LapCoeffs,ExParams);%((R0,th,k);

   if strcmpi(BType,'Chebyshev')
            Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f,dfdn,ChebyshevRange);
        elseif strcmpi(BType,'Fourier')
            Basis = Tools.Basis.FourierBasis.BasisHelper(f,dfdn);
        end

for n=1:6 %run different grids
tic
    p=3;%3;
	Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
		
    Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
    %     WaveNumberClsHandle = @Tools.WaveNumber.ConstantWaveNumber; %@WaveNumberElliptical;
    %     WaveNumberAddParams = k;% struct('k0',k0,'r0',NHR);

    %     sigma=0;
    %     coeffAfunc = @(x)  1;%0.25*sin(x/100).^2+1.3;
    %     coeffBfunc = @(x)  1;%0.2*sin(x/10).^2+1.3;
    
    CoeffsClsHandle = @Tools.Coeffs.ConstLapCoeffs;
    %  CoeffsClsAddParams = struct('a',coeffAfunc(Grid.X),'b', coeffBfunc(Grid.Y),'sigma',sigma);
    CoeffsClsAddParams = LapCoeffs;
    
    if strcmpi(ScatType,'ellipse')
        ScattererClsHandle  = @Tools.Scatterer.EllipticScatterer;%???EllipticScatterer
        ScattererAddParams  = struct('Eta0',Eta0,'FocalDistance',FocalDist);
    elseif strcmpi(ScatType,'circle')
        ScattererClsHandle  = @Tools.Scatterer.PolarScatterer;
        ScattererAddParams  = struct('r0',R0,'ExpansionType',33);
    elseif strcmpi(ScatType,'submarine')
        ScattererClsHandle  = @Tools.Scatterer.SubmarineScatterer;
        ScattererAddParams  = struct('ellipse',ellipse,'tower',tower,'ExpansionType',25);
    end

    IntPrb =  Solvers.InteriorLaplacianHomoSolver ... 
        (Basis,Grid,CoeffsClsHandle,CoeffsClsAddParams,ScattererClsHandle,ScattererAddParams);
    
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
        xiex = Exact(IntPrb.Scatterer.th,LapCoeffs,ExParams2);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
    elseif strcmpi(ScatType,'submarine')
        %ExParams  = struct('ScattererType','submarine','ellipse',ellipse,'tower',tower);
%         xiex = Exact(IntPrb.Scatterer.th,k,ExParams);
        ExParams2 =ExParams;
        ExParams2.r = IntPrb.Scatterer.r;
        xiex = Exact(IntPrb.Scatterer.th,k,ExParams2);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);

    end
    
    %xiex = Exact(Eta(Split.GridGamma),phi,k0);
    ebinf(n) =norm(xiex -xi(IntPrb.GridGamma),inf);
    
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
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),LapCoeffs,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
    elseif strcmpi(ScatType,'submarine')
        %ExParams  = struct('ScattererType','submarine','ellipse',ellipse,'tower',tower);
%         xiex = Exact(IntPrb.Scatterer.nrml_th,k,ExParams);

%         ExParams3 =ExParams;
%         ExParams3.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
%         exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);

        ExParams3 =ExParams;
        ExParams3.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);

    end
   
    
    
    t2=toc;
    
    etinf(n) =norm((exact(:)-u(:)),inf);
    fprintf('coeffs=%d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',0,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
    
    
   % figure, plot(x,abs(u(fix(Ny/2),:)),'r+-',x,abs(exact(fix(Ny/2),:)),'bo');
%     saveas(gcf,['inhomo_k=' num2str(k) 'N=' num2str(Nx) '.jpg'],'jpg');
    % plot(x,abs(u(:,fix(Ny/2))),'r+-',x,abs(exact(:,fix(Ny/2))),'bo');

end

Linf=log2(etinf(1:end-1)./etinf(2:end))
Lbinf=log2(ebinf(1:end-1)./ebinf(2:end))

end



function g = Exact(th,coeffs,Params)%(R,th,k)
if strcmpi(Params.ScattererType,'ellipse')
    x = Params.FocalDistance * cosh(Params.eta) .* cos(th);
    %          y = Params.FocalDistance * sinh(Params.eta) .* sin(th);
elseif strcmpi(Params.ScattererType,'circle')
    x = Params.r.* cos(th);
    y = Params.r .* sin(th);
elseif strcmpi(Params.ScattererType,'submarine')
    e = (cos(th).^2/Params.ellipse.a^2 + sin(th).^2/Params.ellipse.b^2);
    r = sqrt((1 + Params.tower.c*sin(th).^Params.tower.p) ./ e);
    
    x = r .* cos(th);
    %         y = r .* sin(th);
    
end

%g = exp(1i.*k.*R.*cos(th));
%g = exp(1i.*k.*x);

% g  = exp(1i.*coeffs.b.*th).*Params.r.^2./coeffs.a;

g = sin(coeffs.a .* x).*sin(coeffs.b .* y);

end


function dne = drExact(th,coeffs,Params)%(R,th,k)
if strcmpi(Params.ScattererType,'ellipse')
    dx = Params.FocalDistance * sinh(Params.eta) .* cos(th);
    %         dy = Params.FocalDistance * cosh(Params.eta) .* sin(th);
elseif strcmpi(Params.ScattererType,'circle')
    dxdr = cos(th);
    dydr = sin(th);
    
    x = Params.r.* cos(th);
    y = Params.r .* sin(th);
    
    
elseif strcmpi(Params.ScattererType,'submarine')
    e = (cos(th).^2/Params.ellipse.a^2 + sin(th).^2/Params.ellipse.b^2);
    r = sqrt((1 + Params.tower.c*sin(th).^Params.tower.p) ./ e);
    %         dr = Params.tower.c * Params.tower.p * cos(phi).*sin(phi).^(Params.tower.p-1)./e ...
    %             - 2*(-1/Params.ellipse.a^2 + 1/Params.ellipse.b^2).*cos(phi).* sin(phi).*r.^2 ./e;
    %         dr = dr/2/r;
    de = (-1/Params.ellipse.a^2 + 1/Params.ellipse.b^2).* sin(2*th);
    dg = Params.tower.c*Params.tower.p*cos(th).*sin(th).^(-1 + Params.tower.p);
    dr = (dg - de.*r.^2)/2./r./e;
    
    dx = cos(th).*dr;
    %         dy = sin(th).*dr;
    
end

%dne  = 2*exp(1i.*coeffs.b.*th).*Params.r./coeffs.a;

dne = coeffs.b .* cos(coeffs.b .*y) .* dydr.* sin(coeffs.a .*x) + coeffs.a.*dxdr.* cos(coeffs.a .*x).* sin(coeffs.b .*y);

end
