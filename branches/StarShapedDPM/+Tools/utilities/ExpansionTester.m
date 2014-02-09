err=[];

%k=3;
ScatType = 'submarine'; 

kolobok=0;
ellipse = 0;
circle=0;

shape='submarine';
if kolobok
    shape='kolobok';
	AR=1.2;
	c=1.2;

	x1=-3;xn=3;
	y1=-3;yn=3;%=-0.7;yn=0.7;%	
elseif ellipse
    shape = 'ellipse';
 	AR=2;
	c=0;

	x1=-2;xn=2;
	y1=-2;yn=2;
   
elseif circle 
    shape = 'circle';
    AR=1;
	c=0;

	x1=-2;xn=2;
	y1=-2;yn=2;

else %submarine
	AR=2;
	c=4;
	%x1=-2.2; xn=2.2;
	%y1=-2.2; yn=2.2;
    
    x1=-3; xn=3;
	y1=-3; yn=3;
end

a=1.8;
b=a/AR;
p=20;
ellipse = struct('a',a,'b',b);
%tower = struct('c',1.2,'p',p);
tower = struct('c',c,'p',p);

fprintf('%s,params: a=%d, b=%d (AR=%d), c=%d, p=%d, %d<x,y<%d \n',shape,a,b,AR,c,p,x1,xn);

NoSource = Tools.Source.SuperHelmholtzSource();  %construct empty source           
ScattererAddParams  = struct('ellipse',ellipse,'tower',tower,'ExpansionType',25);
for k=1;%[0.1,0.5,1,2,3,5,10]
for n=1:5 %run different grids    
    p=5;%3;
	Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
	%dx=Lx/(Nx-1); 	dy=Ly/(Ny-1);
    	
    Grid        = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
    Scatterer   = Tools.Scatterer.SubmarineScatterer(Grid,ScattererAddParams);
    Submarine   = Scatterer.Submarine;
    t           = Scatterer.BasisArg;
        
    WaveNumberClsHandle = @Tools.WaveNumber.ConstantWaveNumber; %@WaveNumberElliptical;
    WaveNumberAddParams.k = k;

    Coeffs = WaveNumberClsHandle(Scatterer.TheScatterer,WaveNumberAddParams);
       
    %initialization of empty class
    Xi0 = Tools.Basis.BasisFunctionWD(); 
    Xi1 = Tools.Basis.BasisFunctionWD();  
   
    [r,rt,rtt,r3t,r4t] = Submarine.Derivatives(t);
    
    x   = r .*cos(t);
    y   = r .*sin(t);       
    
    xt  = -y    + rt.*cos(t);
    yt  =  x    + rt.*sin(t);
    xtt = -2*yt + x + rtt.*cos(t);
    ytt =  2*xt + y + rtt.*sin(t);
    x3t = -3*ytt+ 3*xt + y + r3t.*cos(t);
    y3t =  3*xtt+ 3*yt - x + r3t.*sin(t);
    x4t = -4*y3t+ 6*xtt +4*yt - x + r4t.*cos(t);
    y4t =  4*x3t+ 6*ytt - 4*xt - y + r4t.*sin(t);
    %x5t = -5*y4t+ 10*x3t + 10*ytt - 5*xt - y + r5t.*cos(t);
    %     y5t =  5*x4t+ 10*y3t - 10*xtt - 5*yt + x + r5t.*sin(t);
    %     x6t = -6*y5t+ 15*x4t + 20*y3t - 15*xtt - 6*yt + r6t.*cos(t);
    
    
    u   = exp(1i*k*x);
    ut  = 1i*k*u.*xt;
    utt = 1i*k*(ut.*xt + u.*xtt);
    u3t = 1i*k*(utt.*xt + 2*ut.*xtt + u.*x3t);
    u4t = 1i*k*(u3t.*xt + 3*utt.*xtt + 3*ut.*x3t + u.*x4t);
    %    u5t = 1i*k*(u4t.*xt + 4*u3t.*xtt + 6*utt.*x3t + 4*ut.*x4t + u.*x5t);
    %     u6t = 1i*k*(u5t.*xt + 5*u4t.*xtt + 10*u3t.*x3t + 10*utt.*x4t + 5*ut.*x5t + u.*x6t);
    
    Xi0.xi0        = u;
    Xi0.xi0t       = ut;
    Xi0.xi0tt      = utt;
    Xi0.xi0ttt     = u3t;
    Xi0.xi0tttt    = u4t;
    %   Xi0.xi0tttttt  = u6t;

       hs      = sqrt(r.^2+rt.^2);
    hst     = (rt.*rtt + r.*rt)./hs;
    hstt    = (rt.^2 + rtt.*r + rtt.^2 + rt.*r3t - hst.^2)./hs;
    hs3t    = (3*rt.*rtt + r3t.*r + 3*rtt.*r3t + rt.*r4t - 3*hst.*hstt)./hs;
    
    ux = 1i*k*u;
    n_x = yt./hs; %(rt.*sin(t) + x)./hs;
    
    n_xt    = (ytt - n_x.*hst)./hs;
    n_xtt   = (y3t - 2*n_xt.*hst - n_x.*hstt)./hs;
    n_x3t   = (y4t - 3.*n_xtt.*hst - 3.*n_xt.*hstt - n_x.*hs3t )./hs;
    
    uxt     = 1i*k*ut;
    uxtt    = 1i*k*utt;
    ux3t    = 1i*k*u3t;
    ux4t    = 1i*k*u4t;
    
    un      = ux.*n_x; % + uy.*n_y,  uy=0
    unt     = uxt.*n_x + ux.*n_xt;
    untt    = uxtt.*n_x + 2*uxt.*n_xt + ux.*n_xtt;
    un3t    = ux3t.*n_x + 3*uxtt.*n_xt + 3*uxt.*n_xtt + ux.*n_x3t;
    

%     un   = ur;
%     uns  = urt./hs;
%     un2s = (urt - uns.*hst)./hs;
%     un3s = (urtt - 2*un2s.*hst -uns.*hstt)./hs;
    
    Xi1.xi0        = un;
    Xi1.xi0t       = unt;
    Xi1.xi0tt      = untt;
    Xi1.xi0ttt     = un3t;
    %     Xi1.xi0tttt    = ur4t;
    %     Xi1.xi0tttttt  = ur6t;
    
    if 1 %test
        xi = Scatterer.Expansion(Xi0,Xi1,NoSource,Coeffs);
        MinDn(n,1) = min(abs(Scatterer.dn(:)));
		MaxDn(n,1) = max(abs(Scatterer.dn(:)));
    else
        xi = u + Scatterer.dn.*un;% + ((obj.dn.^2)/2).*u_nn + obj.dn.^3./6.*u_nnn + obj.dn.^4./24.*u_nnnn;
    end
             
    % the following x,y are input test locations of grid gamma
    TstX = Scatterer.r.*cos(Scatterer.th);
    TstY = Scatterer.r.*sin(Scatterer.th);
    
    % the following x,y are result of finding the normal to the body from
    % the prefious x,y's, we want to verify if they appears to be the same
    % as in the input...
    dtds    = sqrt(r.^2+rt.^2);
    NrmlX = x + Scatterer.dn.*yt./dtds;
    NrmlY = y - Scatterer.dn.*xt./dtds;
        
    
    
    
    
    exact = exp(1i*k*TstX);
    
	tmp = exact(:) - xi(:);
    err(n,1) = norm( tmp ,inf);
	ErrCords(n).GG_i=find(abs(tmp)==err(n));
    ErrCords(n).dn = Scatterer.dn(ErrCords(n).GG_i);
    ErrCords(n).r = Scatterer.r(ErrCords(n).GG_i);
    ErrCords(n).th = Scatterer.th(ErrCords(n).GG_i);
    ErrCords(n).nrml_th = Scatterer.nrml_th(ErrCords(n).GG_i);
    [ErrCords(n).Grid_i,ErrCords(n).Grid_j] = ind2sub([Nx,Ny],Scatterer.GridGamma(ErrCords(n).GG_i));
    ErrCords(n).diffX = abs(NrmlX(ErrCords(n).GG_i) - TstX(ErrCords(n).GG_i));
    ErrCords(n).diffY = abs(NrmlY(ErrCords(n).GG_i) - TstY(ErrCords(n).GG_i));
    
    N(n,1) = Nx;
	
% 	Scat.r = Scatterer.r;
% 	Scat.t = Scatterer.th;
% 	Scat.nt= Scatterer.nrml_th;
% 	Scat.dn= Scatterer.dn;
% 	
% 	ETvars{n} = Scat;
        
end

% save('ETvars.mat','ETvars');
%    err=err.';  
%    MinDn = MinDn.';
%    N=N.';
    Conv = log2(err(1:end-1)./err(2:end));
%    spase = zeros(size(err));
%    spase(1:end) = ' ';
%    [num2str(N), spase, num2str(MinDn), spase, num2str(MaxDn),spase, num2str(err),spase,num2str([NaN;Conv])]
   
   FormatSpectForNUM2STR =  '%-6.4e';
   FS = FormatSpectForNUM2STR;

   Conv = [NaN;Conv];
   	for j=1:numel(err)
 		fprintf('k=%-3.0d\t N=%-4.0d\t Min |n|=%-5.4d\t Max |n|=%-5.4d\t err=%-10.6e \t rate=%-4.2f err_thOnScat=%s\t \t err_th=%s\t err_r=%s\t err_dist=%s\t |x-xn|=%s\t |y-yn|=%s\t\n',...
                        ... %err loc (%s   ;  %s)\n',...
                        k,N(j),MinDn(j),MaxDn(j),err(j),Conv(j),...
                        num2str(ErrCords(j).nrml_th(1),FS),...
                        num2str(ErrCords(j).th(1),FS),...
                        num2str(ErrCords(j).r(1),FS),...
                        num2str(ErrCords(j).dn(1),FS),...
                        num2str(ErrCords(j).diffX(1),FS),...
                        num2str(ErrCords(j).diffY(1),FS)...
                        ... %num2str( ErrCords(j).Grid_i'),...
                        ... %num2str(ErrCords(j).Grid_j')...
                        );
    end
    fprintf('\n');
   
end   
   
   
   

