function RunInterior
	%inhomogenious helmholtz equation with elliptical scatterer and variable wavenumber
	% designed for elliptical scatterer and wavenumbers only

x1=-1.2;xn=1.2;
y1=-0.7;yn=0.7;
Lx=xn-x1;Ly=yn-y1;
ebinf=[];etinf=[];

NHR = 1.6;%1.6;
% k=1;
a=1;
b=0.5;
%for b=0.1:0.1:0.9
    fprintf('Internal Inhomogeneous problem inside ellipse of AR=%d \n',a/b);
FocalDist = sqrt(a^2-b^2);
Eta0 =acosh(1/FocalDist);

for k0 = 1% [1,3,5] %[1,5,10,15,20,25]
    
    f   =@(phi) Exact(FocalDist,Eta0,phi,k0,NHR);
    dfdn=@(phi)  detaExact(FocalDist,Eta0,phi,k0,NHR);

Basis = Tools.Basis.FourierBasis.BasisHelper(f,dfdn);
  %  [cn0,cn1ex,M] = FourierCoeff(f,dfdn);
    % [cn0,cn1ex] = FourierCoeff2(f,dfdr,M);

for n=1:4 %run different grids
tic
	%build grid
% 		Nx=2.^(n+1)+5;	Ny=2.^(n+1)+5;
p=3;%3;
	Nx=2.^(n+p);	Ny=2.^(n+p);
	%dx=Lx/(Nx-1); 	dy=Ly/(Ny-1);
    
	
    %BasisIndices        = -M:M;
    
    Grid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
    WaveNumberClsHandle = @Tools.WaveNumber.WaveNumberElliptical;
    WaveNumberAddParams = struct('k',k0,'r0',NHR);
    ScattererClsHandle  = @Tools.Scatterer.EllipticScatterer;%Internal
    ScattererAddParams  = struct('Eta0',Eta0,'FocalDistance',FocalDist);
    Source              = @Tools.Source.HelmholtzSource;
    
    
%     (BasisIndices,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererC
%     lsHandle,ScattererAddParams,Source)
    
    IntPrb = Solvers.InteriorSolver ... 
        (Basis,Grid,WaveNumberClsHandle,WaveNumberAddParams,ScattererClsHandle,ScattererAddParams,Source);%(k0,x1,xn,y1,yn,Nx,Ny,FocalDist,Eta0,M,NHR);
    
    Q0 = IntPrb.Q0;%(:,1:2*M+1);
    Q1 = IntPrb.Q1;%(:,2*M+2:4*M+2);
    Qf = IntPrb.Qf;%(:,2*M+2:4*M+2);
    
    TrGF = IntPrb.TrGF;
%     clear GLW
    
    cn1 =( Q1 \ ( -Q0*Basis.cn0 - TrGF - Qf)) ;
    
        
    xi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
    xi(IntPrb.GridGamma) = ...
        IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1 + IntPrb.Wf(IntPrb.GridGamma);
ebinf(n)=0;
%       xi(Split.GridGamma)=getIxfaceCond(dr,th,k0,IntR);%test
    xiex = Exact(FocalDist,IntPrb.Scatterer.eta,IntPrb.Scatterer.phi,k0,NHR);
    %xiex = Exact(Eta(Split.GridGamma),phi,k0);
    ebinf(n) =norm(xiex -xi(IntPrb.GridGamma),inf);
    
 %xi(Split.GridGamma)=xiex;
    
%%%%%%%%%%%%
%         u = spalloc(Nx,Ny,length(IntPrb.Np));
%         
%         GLW = IntPrb.Solve(xi);
%         
%         u(IntPrb.Np)=xi(IntPrb.Np) - GLW(IntPrb.Np).';
% 
%         u=u + IntPrb.GF;

        u = IntPrb.P_Omega(xi);


%     tmp = A*xi(:);
%     rhs=zeros(cols,rows);
%     rhs(Split.Mp)=tmp(Split.Mp);
%     GLW = A\rhs(:);
%     tmp=xi(:)-GLW;
%     tmp = reshape(tmp, cols, rows);
%     u(Split.Np)=tmp(Split.Np);
%     u=u + GF;

   % ebinf(n) =norm(xi(Split.GridGamma) -u(Split.GridGamma),inf);
    
    %%%%%%%%%%%%%%
    
   t1=toc; 

   % % % % % % % % % % % % % % % %
    % Comparison
    % % % % % % % % % % % % % % % %
    
    
    
    
    
    exact = zeros(Nx,Ny);
    exact(IntPrb.Scatterer.Np) = Exact(FocalDist,IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np),k0,NHR);
       
    t2=toc;
    
    etinf(n) =norm(exact(:)-u(:),inf);
    fprintf('k=%d,M=%d,N=%-10dx%-10d\t ebinf=%d\tetinf=%d\ttimeA=%d\ttimeE=%d\n',k0,Basis.M, Nx,Ny,full(ebinf(n)),full(etinf(n)),t1,t2-t1);
    
    
   % figure, plot(x,abs(u(fix(Ny/2),:)),'r+-',x,abs(exact(fix(Ny/2),:)),'bo');
%     saveas(gcf,['inhomo_k=' num2str(k) 'N=' num2str(Nx) '.jpg'],'jpg');
    % plot(x,abs(u(:,fix(Ny/2))),'r+-',x,abs(exact(:,fix(Ny/2))),'bo');

end

Linf=log2(etinf(1:end-1)./etinf(2:end))
Lbinf=log2(ebinf(1:end-1)./ebinf(2:end))
end
%end
end


function e = Exact(FocDist,eta,phi,k0,r0)

    k = Tools.WaveNumber.WaveNumberElliptical.kkn(FocDist,eta,phi,k0,r0);

    x = FocDist*cosh(eta).*cos(phi);
    e = exp(1i.*k.*x); 
end

function dne = detaExact(FocDist,eta,phi,k0,r0)

    [k,kn] = Tools.WaveNumber.WaveNumberElliptical.kkn(FocDist,eta,phi,k0,r0);
    
    e = Exact(FocDist,eta,phi,k0,r0);
    dne = 1i.*FocDist.*cos(phi).*(kn.*cosh(eta)+k.*sinh(eta)).*e;
end
