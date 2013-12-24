% use when stop on brake point in scatterer

%% ellipse mode

 r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
  mesh(X,Y,ones(size(X)));

tPhi=0:0.001:2*pi;
hold
plot(obj.FocalDistance*cosh(obj.Eta0).*cos(tPhi),obj.FocalDistance*sinh(obj.Eta0).*sin(tPhi))


plot(obj.FocalDistance*cosh(obj.eta).*cos(obj.phi),obj.FocalDistance*sinh(obj.eta).*sin(obj.phi),'b.')

  
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
% saveas(gcf,'GridGammaEllipse2.jpg','jpg')


%% circle mode
polar=0;
if polar == 1
r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
 
else
    X=obj.Grid.X();
    Y=obj.Grid.Y();
end

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
 view(0,-90)

% use after the figure looks fine
%saveas(gcf,'GridGammaCrcl2.jpg','jpg')

%% submarine mode
polar=0;
if polar == 1
 r = obj.Grid.r;
 th = [obj.Grid.theta,2*pi];
 [R,Th]=meshgrid(r,th);

 X=R.*cos(Th);
 Y=R.*sin(Th);
 
else
    X=obj.Grid.X();
    Y=obj.Grid.Y();
end

  mesh(X,Y,ones(size(X)));

hold
tPhi=0:0.001:2*pi;
s=obj.Submarine.Derivatives(tPhi);

plot(s.*cos(tPhi),s.*sin(tPhi))


plot(obj.r.*cos(obj.th),obj.r.*sin(obj.th),'b.')

tPhi=obj.nrml_th;%0:0.001:2*pi;
s=obj.Submarine.Derivatives(tPhi);

  
h=gca;
                set(h,'Color','none');
                set(h,'Visible','off');


hold off

for m=1:25
    line([s(m).*cos(obj.nrml_th(m)),obj.r(m).*cos(obj.th(m))], ...
        [s(m).*sin(obj.nrml_th(m)),obj.r(m).*sin(obj.th(m))]);
%      line([s(m).*cos(obj.nrml_th(m)),0], ...
%         [s(m).*sin(obj.nrml_th(m)),0],'Color','r');
end

%axis([-2,2,-2.2,2.2])
axis equal
 view(0,-90)
% use after the figure looks fine
% saveas(gcf,'GridGammaEllipse2.jpg','jpg')
