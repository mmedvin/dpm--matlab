
x1=-3;xn=3;
y1=-3;yn=3;%=-0.7;yn=0.7;%
Lx=xn-x1;Ly=yn-y1;

k=1;
ScatType = 'submarine'; 

% AR=2;

kolobok=0;
ellipse = 0;
circle=0;

if kolobok
	AR=1.2;
	c=1.2;

	x1=-3;xn=3;
	y1=-3;yn=3;%=-0.7;yn=0.7;%	
elseif ellipse
 	AR=1.2;
	c=0;

	x1=-2;xn=2;
	y1=-2;yn=2;
   
elseif circle 
    AR=1;
	c=0;

	x1=-2;xn=2;
	y1=-2;yn=2;

else %submarine
	AR=2;
	c=4;
	x1=-2.2; xn=2.2;
	y1=-2.2; yn=2.2;
end

a=1.8;
b=a/AR;
ellipse = struct('a',a,'b',b);
%tower = struct('c',1.2,'p',20);
tower = struct('c',c,'p',20);



NoSource = Tools.Source.SuperHelmholtzSource();  %construct empty source           
ScattererAddParams  = struct('ellipse',ellipse,'tower',tower,'ExpansionType',25);

    n=1;
    p=5;%3;
	Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
	%dx=Lx/(Nx-1); 	dy=Ly/(Ny-1);
    	
    Grid        = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
    
    
           
    TestParams.dn = 1./[1,2,4,8,16,32,64,128,256,512];

    for n = 1:2:10
        TestParams.Angle = pi/n;
    
    Scatterer   = Tools.Scatterer.TesterScatterer(Grid,ScattererAddParams,TestParams); %Tools.Scatterer.SubmarineScatterer(Grid,ScattererAddParams);
    Submarine   = Scatterer.Submarine;
    t           = Scatterer.BasisArg;
    
        
    WaveNumberClsHandle = @Tools.WaveNumber.ConstantWaveNumber; %@WaveNumberElliptical;
    WaveNumberAddParams = k;

    Coeffs = WaveNumberClsHandle(Scatterer.TheScatterer,WaveNumberAddParams);
       
    %initialization of empty class
    Xi0 = Tools.Basis.BasisFunctionWD(); 
    Xi1 = Tools.Basis.BasisFunctionWD();  
   
    [r,rt,rtt,r3t,r4t,r5t] = Submarine.Derivatives(t);
    
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
    x5t = -5*y4t+ 10*x3t + 10*ytt - 5*xt - y + r5t.*cos(t);
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

    %     xr = cos(t);
    %     xrt = -sin(t);
    %     xrtt = -xr;
    %     xr3t = -xrt;
    %     xr4t = - xrtt;
    %
    %     xr  = rt.*cos(t);
    %
    %     xrt  = -rt.*sin(t) + rtt.*cos(t)  ;
    %     xrtt = -xr - 2*rtt.*sin(t) + r3t.*cos(t) ;
    %     xr3t = -xrt - 2*cos (t).*rtt - 3*r3t.*sin (t) + r4t.*cos (t);
    %     xr4t = -2*xrtt - xr - 4*cos (t).*r3t - 4*sin (t).*r4t + cos (t).*r5t;
    
    %     xr = xt./rt;
    %     xrt = (xtt - xr.*rtt)./rt;
    %     xrtt = (x3t - 2*xrt.*rtt - xr.*r3t)./rt;
    %     xr3t = (x4t - 3*xrtt.*rtt - 3*xrt.*r3t - xr.*r4t)./rt;
    %     xr4t = (x5t - 4*xr3t.*rtt - 6*xrtt.*r3t - 4*xrt.*r4t - xr.*r5t)./rt;
    
    
    %     dudr = dudx*dxdt*dtdr
    %     ur = 1i*k*u .* xt ./rt; ==> xr = xt./rt
    
    %     xr = xt./rt;
    %     xrt = (xtt - xr.*rtt)./rt;
    %     xrtt = (x3t - 2*xrt.*rtt - xr.*r3t)./rt;
    %     xr3t = (x4t - 3*xrtt.*rtt - 3*xrt.*r3t - xr.*r4t)./rt;
    %     xr4t = (x5t - 4*xr3t.*rtt - 6*xrtt.*r3t - 4*xrt.*r4t - xr.*r5t)./rt;
    %
    %
    %     ur      = 1i*k*u.*xr;
    %     urt     = 1i*k*(ut.*xr + u.*xrt);
    %     urtt    = 1i*k*(utt.*xr + 2*ut.*xrt + u.*xrtt);
    %     ur3t    = 1i*k*(u3t.*xr + 3*utt.*xrt + 3*ut.*xrtt + u.*xr3t);
    %     ur4t    = 1i*k*(u4t.*xr + 4*u3t.*xrt + 6*utt.*xrtt + 4*ut.*xr3t + u.*xr4t);
    %     ur5t    = 1i*k*(u5t.*xr + 5*u4t.*xrt + 10*u3t.*xrtt + 10*utt.*xr3t + 5*ut.*xr4t + u.*xr5t);
    %     ur6t    = 1i*k*(u6t.*xr + 6*u5t.*xrt + 15*u4t.*xrtt + 20*u3t.*xr3t + 15*utt.*xr4t + 6*ut.*xr5t + u.*xr6t);
    
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
    
  if 0  
%     un      = ur./hr;
%     unt     = (urt - un.*hrt)./hr;
%     untt    = (urtt - 2*unt.*hrtt -un.*hrtt)./hr;
%     un3t    = (ur3t - 3*untt.*hrt - 3*unt.*hrtt - un.* hr3t)./hr;

%     x_n = xt./hr;
%     xnt = (xtt - x_n.*hrt)./hr;
%     xntt = (x3t - 2*xnt.*hrt - x_n.*hrtt)./hr;
%     xn3t = (x4t - 3*xntt.*hrt - 3*xnt.*hrtt - x_n.*hr3t)./hr;
% 
%     un      = 1i*k*u.*x_n;
%     unt     = 1i*k*(ut.*x_n + u.*xnt);
%     untt    = 1i*k*(utt.*x_n + 2*ut.*xnt + u.*xntt);
%     un3t    = 1i*k*(u3t.*x_n + 3*utt.*xnt + 3*ut.*xntt + u.*xn3t);    
  end
  
    Xi1.xi0        = un;
    Xi1.xi0t       = unt;
    Xi1.xi0tt      = untt;
    Xi1.xi0ttt     = un3t;
%     Xi1.xi0tttt    = ur4t;
    %     Xi1.xi0tttttt  = ur6t;
        
    
    if 1 %test
        xi = Scatterer.Expansion(Xi0,Xi1,NoSource,Coeffs);
    else
        xi = u + Scatterer.dn.*ur;% + ((obj.dn.^2)/2).*u_nn + obj.dn.^3./6.*u_nnn + obj.dn.^4./24.*u_nnnn;
    end
        
    
    
    TstX = Scatterer.r.*cos(Scatterer.th);
        
    exact = exp(1i*k*TstX);
    
    err = abs(exact(:) - xi(:));
%    indx = find(isnan(tmp));
 %   err(indx) = exact(indx);
    
   % err(n) = norm(tmp,inf);        

   %err=err.';
   Conv = log2(err(1:end-1)./err(2:end));
   
  %  dn = 1./[1,2,4,8,16,32,64,128,256,512,1024,2048,4096].';
   
   
   
   spase = zeros(size(err));
   spase(1:end) = ' ';
   Conv = [NaN;Conv];
   for j=1:numel(TestParams.dn)
       %fprintf('%-13.4e \t %s \t %-13.4e \n',dn(j),err(j),Conv(j));
       fprintf('Angle=pi/%-6.0d\tn=1/%-6.0d \t err=%-12.6e \t rate=%-10.2f \n',pi/TestParams.Angle,1/TestParams.dn(j),err(j),Conv(j));
   end
   %fprintf('n=%-13.4e \t err=%-15.6e \t rate=%-13.4e \n',dn,err,[NaN;Conv]);
   %[num2str(dn),spase, num2str(err),spase,num2str([NaN;Conv])]
   fprintf('\n');
    end