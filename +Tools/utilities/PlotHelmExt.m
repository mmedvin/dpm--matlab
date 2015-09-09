%%

R  = ones(size(Grid.R)).*NaN;
Th = ones(size(Grid.R)).*NaN;
Nm = ExtPrb.Scatterer.Nm;
R(Nm) = Grid.R(Nm);
Th(Nm)= Grid.Theta(Nm);

R(:,Nth)=R(:,Nth-1);
Th(:,Nth)=2*pi;

XExt = (R .* cos(Th));
YExt = (R .* sin(Th));

tExtu= ones(size(Grid.R)).*NaN;
tExtu(Nm)=u(Nm);
tExtu(:,Nth)=tExtu(:,1);


 th=0:0.0001:2*pi;
                
 figure(1)
 pcolor(XExt,YExt,abs(full(tExtu)))
 
 colormap jet
 title('scattered field, abs');
 axis equal
 axis off;
 view(2);
 shading flat;
 grid off;
 h=gca;
 set(h,'Color','none');
 set(h,'Visible','off');
 colorbar
 
 hold on;
 %plot(a*cos(th),b*sin(th),'k','LineWidth',2)
 plot(Parameterization.XHandle.Derivatives(th),Parameterization.YHandle.Derivatives(th),'k','LineWidth',2)
 
 hold off;
 
 figure(2)
  pcolor(XExt,YExt,real(full(tExtu)))
 
 colormap jet
 title('scattered field, real');
 axis equal
 axis off;
 view(2);
 shading flat;
 grid off;
 h=gca;
 set(h,'Color','none');
 set(h,'Visible','off');
 colorbar
 
 hold on;
 plot(Parameterization.XHandle.Derivatives(th),Parameterization.YHandle.Derivatives(th),'k','LineWidth',2)
 hold off;
 
 figure(3)
  pcolor(XExt,YExt,imag(full(tExtu)))
 
 colormap jet
 title('scattered field, imag');
 axis equal
 axis off;
 view(2);
 shading flat;
 grid off;
 h=gca;
 set(h,'Color','none');
 set(h,'Visible','off');
 colorbar
 
 hold on;
 plot(Parameterization.XHandle.Derivatives(th),Parameterization.YHandle.Derivatives(th),'k','LineWidth',2)
 hold off;
 
 % tmp=abs([tIntu,tExtu]);
 % fprintf('min=%d \t max=%d', full(min(tmp(:))), full(max(tmp(:))));
 %
 % saveas(figure(1),[filename 'abs.jpg'],'jpg')
 % saveas(figure(2),[filename 'real.jpg'],'jpg')
 % saveas(figure(3),[filename 'imag.jpg'],'jpg')
 
 %fprintf('press key');
 %soundsc(sin(2000*(0:0.001:2*pi)).^24);
 %pause
 %fprintf('.');