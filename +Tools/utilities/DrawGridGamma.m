% use when stop on brake point in scatterer

%% starshaped body mode
figure
hold on

try
    X=obj.Grid.X;
    Y=obj.Grid.Y ;
catch

 r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
 
end

 mesh(X,Y,ones(size(X)), 'EdgeAlpha',0.7);%'LineStyle',':');%,'EdgeAlpha',0.7
%  light('FaceColor', 'interp')

	tPhi=0:0.001:2*pi;
	plot(obj.XHandle.Derivatives(tPhi),obj.YHandle.Derivatives(tPhi),'LineWidth',2)
	
	plot(obj.r.*cos(obj.th),obj.r.*sin(obj.th),'b.','MarkerSize',10)
	
    
	
	axis equal
	
	h=gca;
	set(h,'Color','none');
	set(h,'Visible','off');
    %shading flat;
    
	%view(0,-90) %kite
    view(180,-90)
	
	hold off
% 	
% 	for m=1:numel(obj.GridGamma)
% 		line([obj.XHandle.Derivatives(obj.nrml_t(m)),obj.r(m).*cos(obj.th(m))], ...
% 					[obj.YHandle.Derivatives(obj.nrml_t(m)), obj.r(m).*sin(obj.th(m))]);
% 	end

%%
%%%%%%%%%%%%% draw wave number 


Scatt = struct('r',obj.Grid.R);
WN = Tools.Coeffs.WaveNumberPolarR(Scatt,WaveNumberAddParams);

k = WN.k;
kmp=NaN*zeros(size(k));
%kmm=NaN*zeros(size(k));
kmp(obj.Scatterer.Mp) = k(obj.Scatterer.Mp);
%kmp(obj.Scatterer.Mm) = k(obj.Scatterer.Mm);

try
    X=obj.Grid.X;
    Y=obj.Grid.Y ;
catch

 r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
	
	end

XMp = NaN*zeros(size(k));
YMp = NaN*zeros(size(k));
XMm = NaN*zeros(size(k));
YMm = NaN*zeros(size(k));

XMp(obj.Scatterer.Mp) = X(obj.Scatterer.Mp);
YMp(obj.Scatterer.Mp) = Y(obj.Scatterer.Mp);

XMm(obj.Scatterer.Mm) = X(obj.Scatterer.Mm);
YMm(obj.Scatterer.Mm) = Y(obj.Scatterer.Mm);


figure 
%mesh(X,Y,k)
%colormap default, 
colormap jet
m=mesh(XMp,YMp,kmp,'FaceColor','interp'); hold on, 
%m=mesh(XMm,YMm,kmm,'FaceColor','flat');
m=mesh(XMm,YMm,1*ones(size(X)),'FaceColor','flat');
%hold off, camlight, 
lighting gouraud, 
%xlabel('x'), ylabel('y'),  
%set(m,'facecolor','cyan','edgecolor','none'),
%set(m,'edgecolor','none'),
%alpha(0.5)

%axis equal
	    
h=gca;
set(h,'Color','none');
%set(h,'Visible','off');

%%%%%%%%%%%%%%%

%% ellipse mode
figure
try
    X=obj.Grid.X;
    Y=obj.Grid.Y ;
catch

 r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
  mesh(X,Y,ones(size(X)));
end
tPhi=0:0.001:2*pi;
hold
plot(obj.FocalDistance*cosh(obj.Eta0).*cos(tPhi),obj.FocalDistance*sinh(obj.Eta0).*sin(tPhi))


plot(obj.FocalDistance*cosh(obj.eta).*cos(obj.phi),obj.FocalDistance*sinh(obj.eta).*sin(obj.phi),'b.')

  
h=gca;
                set(h,'Color','none');
                set(h,'Visible','off');


hold off

for m=1:numel(obj.GridGamma)
try
         line([obj.XHandle.Derivatives(obj.nrml_t(m)),obj.FocalDistance*cosh(obj.eta(m)).*cos(obj.phi(m))], ...
         [obj.YHandle.Derivatives(obj.nrml_t(m)),obj.FocalDistance*sinh(obj.eta(m)).*sin(obj.phi(m))]);
catch
    line([obj.FocalDistance*cosh(obj.Eta0).*cos(obj.phi(m)),obj.FocalDistance*cosh(obj.eta(m)).*cos(obj.phi(m))], ...
         [obj.FocalDistance*sinh(obj.Eta0).*sin(obj.phi(m)),obj.FocalDistance*sinh(obj.eta(m)).*sin(obj.phi(m))]);
end
end


axis equal

% use after the figure looks fine
% saveas(gcf,'GridGammaEllipse2.jpg','jpg')


%% circle mode
r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
  mesh(X,Y,ones(size(X)));

tPhi=0:0.001:2*pi;
hold
plot(obj.r0.*cos(tPhi),obj.r0.*sin(tPhi))


plot(obj.r.*cos(obj.th),obj.r.*sin(obj.th),'b.')

  
h=gca;
                set(h,'Color','none');
                set(h,'Visible','off');


hold off

% for m=1:25
%     line([obj.FocalDistance*cosh(obj.Eta0).*cos(obj.phi(m)),obj.FocalDistance*cosh(obj.eta(m)).*cos(obj.phi(m))], ...
%         [obj.FocalDistance*sinh(obj.Eta0).*sin(obj.phi(m)),obj.FocalDistance*sinh(obj.eta(m)).*sin(obj.phi(m))]);
% end


axis equal

% use after the figure looks fine
%saveas(gcf,'GridGammaCrcl2.jpg','jpg')

%% submarine mode

 r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
  mesh(X,Y,ones(size(X)));
hold
tPhi=0:0.001:2*pi;
s=obj.Submarine.Derivatives(tPhi);

plot(s.*cos(tPhi),s.*sin(tPhi))


plot(obj.r.*cos(obj.th),obj.r.*sin(obj.th),'b.')

  
h=gca;
                set(h,'Color','none');
                set(h,'Visible','off');


hold off

% for m=1:25
%     line([obj.FocalDistance*cosh(obj.Eta0).*cos(obj.phi(m)),obj.FocalDistance*cosh(obj.eta(m)).*cos(obj.phi(m))], ...
%         [obj.FocalDistance*sinh(obj.Eta0).*sin(obj.phi(m)),obj.FocalDistance*sinh(obj.eta(m)).*sin(obj.phi(m))]);
% end


axis equal
 view(0,-90)
% use after the figure looks fine
% saveas(gcf,'GridGammaEllipse2.jpg','jpg')
