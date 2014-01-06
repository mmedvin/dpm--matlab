% this one tests the finding closest distance...
% here we use a data saved in a run of ExpansionTester ETvars
% ETvars is a cell of structures, each structure consist arrays of 
% (r,t) - are coordinates, aka (r,theta), of a point on a grid gamma
% nt is a angle in a normal direction from (r,t) to the body 
% dn is a distance from  (r,t) to (R(nt), nt) where R describe the body curve

load('ETvars')

x1=-3;xn=3;
y1=-3;yn=3;%=-0.7;yn=0.7;%
Lx=xn-x1;Ly=yn-y1;

k=1;
ScatType = 'submarine'; 

kolobok=0;
ellipse = 0;
circle=1;

if kolobok
	AR=1.2;
	c=1.2;

	x1=-3;xn=3;
	y1=-3;yn=3;%=-0.7;yn=0.7;%	
elseif ellipse
 	AR=2;
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
    
    
           
    %TestParams.dn = 1./[1,2,4,8,16,32,64,128,256,512,1024];

	M = length(ETvars);
	Save = cell(M,1);
	
	for m=1:M
		Tester=ETvars{m};
		L = length(ETvars{m}.nt);
		R = zeros(L,1);
		T=R;
		for n = 1:L
			TestParams.dn		= Tester.dn(n);
			TestParams.Angle	= Tester.nt(n);
			
			Scatterer   = Tools.Scatterer.TesterSubmarineScatterer(Grid,ScattererAddParams,TestParams); %Tools.Scatterer.SubmarineScatterer(Grid,ScattererAddParams);
			if 0
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
				else
					xi = u + Scatterer.dn.*un;% + ((obj.dn.^2)/2).*u_nn + obj.dn.^3./6.*u_nnn + obj.dn.^4./24.*u_nnnn;
				end
				
				TstX = Scatterer.r.*cos(Scatterer.th);
				
				exact = exp(1i*k*TstX);
				
				err = abs(exact(:) - xi(:));
				
				Conv = log2(err(1:end-1)./err(2:end));
				
				Conv = [NaN;Conv];
				for j=1:numel(TestParams.dn)
					fprintf('Angle=pi/%-6.0d\tn=1/%-6.0d \t err=%-12.6e \t rate=%-10.2f \n',pi/TestParams.Angle,1/TestParams.dn(j),err(j),Conv(j));
				end
			end
			
			R(n) = Scatterer.r;
			T(n) = Scatterer.th;
		end
		
		Scat.r = R;
		Scat.t = T;
		Save{m}=Scat;
		
		%
		% 			figure(1), plot(TesterScat.r - R);
		% 			figure(2), plot(TesterScat.t - T);
		% 			pause;
		%fprintf('\n');
	end
	
	save('Save.mat','Save');