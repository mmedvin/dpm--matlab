
GridN=512;

k           = 5;
IncAng      = pi/5;
BC          = Tools.Enums.BoundaryConditions.Neumann;%Dirichlet;%;

Scat1Type    = Tools.Enums.Scatterer.StarShaped;
Scat2TypeShifted    = Tools.Enums.Scatterer.StarShaped;

HankOrPlane = 'PlaneWave';% 'PlaneWave' or 'Hankel'
HankelIndex = 3; HankelType = 2;

%Ellipse Axes:
a=1;
b=a;

%Ring radiuses:
r0=0.3; r1=2.2; %r0=0.8; r1=2.2;
%r0=0.8; r1=2.2;

BasisType = Tools.Enums.Basis.Fourier;

Shift  =[0.3,0.0];
NoShift=[0.0,0.0];
Parameterization1           = Tools.Parameterizations.ParametricEllipse2(struct('a',a,'b',b,'xcenter',NoShift(1),'ycenter',NoShift(2),'rotation',0));
Parameterization2Shifted    = Tools.Parameterizations.ParametricEllipse2(struct('a',a,'b',b,'xcenter',Shift(1)  ,'ycenter',Shift(2)  ,'rotation',0));

[Scatterer1Handle,Scatterer1Params,Extension1,ExParams1] = Scat1Type.Helper         (a,b,HankelIndex,HankelType,Parameterization1);
[Scatterer2Handle,Scatterer2Params,Extension2,ExParams2] = Scat2TypeShifted.Helper  (a,b,HankelIndex,HankelType,Parameterization2Shifted);


fprintf('ShiftingTest with Scattering Problem about %s\n',Scat1Type);
Scat1Type.Print(a,b,HankelIndex,HankelType,Parameterization1);
fprintf('and about %s\n',Scat2TypeShifted)
Scat2TypeShifted.Print(a,b,HankelIndex,HankelType,Parameterization2Shifted);
fprintf('Polar Grid in Ring between r0=%-5.2f and r1=%-5.2f , GridN = %d\n ',r0,r1,GridN );

fprintf('BC:%s, Basis is %s, k=%f \n',BC.toString,BasisType.toString,k);

Grid = Tools.Grid.PolarGrids(r0,r1,GridN,GridN);
WaveNumberHndl = @Tools.Coeffs.ConstantWaveNumber;
WaveNumberHndlParams = struct('k',k,'r0',1.6);
SetupCommon = struct('Grid', Grid, 'CoeffsHandle', WaveNumberHndl, 'CoeffsParams', WaveNumberHndlParams,'CollectRhs', 0); 

f1      = @(phi) -Uinc    (ExParams1,phi,IncAng,k);
df1dn   = @(phi) -detaUinc(ExParams1,phi,IncAng,k);

Basis1 = BasisType.Helper(f1,df1dn);%,[1e-06,1e-06]);
Setup1 = SetupCommon;
Setup1.Basis = Basis1;
Setup1.ScattererHandle = Scatterer1Handle;
Setup1.ScattererParams = Scatterer1Params;
Setup1.Extension       = Extension1;
Setup1.ExtensionParams =[];

f2      = @(phi) -Uinc    (ExParams2,phi,IncAng,k);
df2dn   = @(phi) -detaUinc(ExParams2,phi,IncAng,k);

Basis2 = BasisType.Helper(f2,df2dn);%,[1e-06,1e-06]);

Setup2Shited = SetupCommon;
Setup2Shited.Basis = Basis2;
Setup2Shited.ScattererHandle = Scatterer2Handle;
Setup2Shited.ScattererParams = Scatterer2Params;
Setup2Shited.Extension       = Extension2;
Setup2Shited.ExtensionParams =[];


[u1,ExtPrb1,u1Cn0,u1Cn1] = Solve(BC,Setup1);
[u2Shifted,ExtPrb2,u2ShiftedCn0,u2ShiftedCn1] = Solve(BC,Setup2Shited);
KincHat = [cos(IncAng-pi),sin(IncAng-pi)];
arg = KincHat(1)*Shift(1) + KincHat(2)*Shift(2); 
u3=u1*exp(-1i*k*arg);
%norm(u3(ExtPrb1.GridGamma)-u2(ExtPrb2.GridGamma),inf)

% figure
% plot(ExtPrb1.Scatterer.th,real(u1(ExtPrb1.GridGamma)),'gv',ExtPrb1.Scatterer.th,real(u3(ExtPrb1.GridGamma)),'r+',ExtPrb2.Scatterer.th,real(u2(ExtPrb2.GridGamma)),'bo')
% legend('u1','u3=u1_shifted','u2shifted')
% figure, 
% plot(ExtPrb1.Scatterer.th,imag(u1(ExtPrb1.GridGamma)),'gv',ExtPrb1.Scatterer.th,imag(u3(ExtPrb1.GridGamma)),'r+',ExtPrb2.Scatterer.th,imag(u2(ExtPrb2.GridGamma)),'bo')
% legend('u1','u3=u1_shifted','u2shifted')

if 0
    warning('do the same test, but using basis instead of interpolation');
    x = cos(ExtPrb1.Scatterer.th).*ExtPrb1.Scatterer.r;
    y = sin(ExtPrb1.Scatterer.th).*ExtPrb1.Scatterer.r;
    tt= linspace(0,2*pi,1000);
    nx = Parameterization1.XHandle.Derivatives(tt);
    ny = Parameterization1.YHandle.Derivatives(tt);
    m3 = griddata(x,y,full(u3(ExtPrb1.GridGamma)),nx,ny);
    
    x = cos(ExtPrb2.Scatterer.th).*ExtPrb2.Scatterer.r;
    y = sin(ExtPrb2.Scatterer.th).*ExtPrb2.Scatterer.r;
    tt= linspace(0,2*pi,1000);
    nx = Parameterization2Shifted.XHandle.Derivatives(tt);
    ny = Parameterization2Shifted.YHandle.Derivatives(tt);
    m2 = griddata(x,y,full(u2(ExtPrb2.GridGamma)),nx,ny);
end
if 1
    tt= linspace(0,2*pi,1000);
    u2OnGamma = Assemble(tt, Basis2.Handle,Basis2.Indices0,u2ShiftedCn0);
    Shiftedu1OnGamma = Assemble(tt, Basis1.Handle,Basis1.Indices0,u1Cn0)*exp(-1i*k*arg);
    
    figure
    sgtitle('u2|_\Gamma and shifted u1|_\Gamma, assembled from Basis');
    subplot(121),  plot(tt,real(Shiftedu1OnGamma),'r+',tt,real(u2OnGamma),'bo'), title('real')
    subplot(122),  plot(tt,imag(Shiftedu1OnGamma),'r+',tt,imag(u2OnGamma),'bo'), title('imag')
    legend('u3=u1_shifted','u2shifted')
    
    
    ttt= linspace(0,2*pi,2000);
    ttt=ttt(1:end-1);
    
    field2.value        = Assemble(ttt, Basis2.Handle,Basis2.Indices0,u2ShiftedCn0);
    field2.normal_deriv = Assemble(ttt, Basis2.Handle,Basis2.Indices1,u2ShiftedCn1);
    
    field3.value        = Assemble(ttt, Basis1.Handle,Basis1.Indices0,u1Cn0)*exp(-1i*k*arg);
    field3.normal_deriv = Assemble(ttt, Basis1.Handle,Basis1.Indices1,u1Cn1)*exp(-1i*k*arg);
    
    [FFPOnGamma_plus2,FFPOnGamma_minus2] = getFFP(field2,setfield(Setup2Shited  ,'theta',ttt),ExParams2);
    [FFPOnGamma_plus3,FFPOnGamma_minus3] = getFFP(field3,setfield(Setup1        ,'theta',ttt),ExParams1);
    
    figure('units', 'normalized', 'position', [0.1 0.1 0.8, 0.4])
    
    phi = FFPOnGamma_minus2.phi_refl;
    xhat = [cos(phi),sin(phi)];
    %IncAng = IncAng-pi;
    KincHat = [cos(IncAng),sin(IncAng)];
    arg = (xhat(1)-KincHat(1))*Shift(1) + (xhat(2)-KincHat(2))*Shift(2); Title = 'FFP minus from u2|_\Gamma and shifted u1|_\Gamma, assembled from Basis';
    
    F2 = FFPOnGamma_minus2.ffp; % computed from the field
    ShitftedF0 = FFPOnGamma_minus3.ffp;
    
    sgtitle(Title);
    subplot(221), plot(phi*180/pi,real(F2) ,'bo',phi*180/pi, real(ShitftedF0),'r+') , title('real' )
    subplot(222), plot(phi*180/pi,imag(F2) ,'bo',phi*180/pi, imag(ShitftedF0),'r+') , title('imag' )
    subplot(223), plot(phi*180/pi,abs(F2)  ,'bo',phi*180/pi, abs(ShitftedF0),'r+')  , title('abs'  )
    subplot(224), plot(phi*180/pi,angle(F2),'bo',phi*180/pi, angle(ShitftedF0),'r+'), title('angle')
    
    diff = F2 - ShitftedF0;
    
    figure('units', 'normalized', 'position', [0.1 0.1 0.8, 0.4])
    sgtitle(Title);
    subplot(221), plot(phi*180/pi,real(diff)), title('real')
    subplot(222), plot(phi*180/pi,imag(diff)), title('imag')
    subplot(223), plot(phi*180/pi,abs(diff)), title('abs')
    subplot(224), plot(phi*180/pi,angle(diff)), title('angle')
    
    fprintf('FFP minus from u on Gamma: max-err(F2-ShiftedF0)=%f, max-err(angle(F2-ShiftedF0))=%f\n', norm(diff,inf)/norm(F2,inf), norm(angle(diff),inf))
    
    
end

if 1
    [UincParams1,Th1] = Scat1Type.UincOnFieldHelper(ExtPrb1,ExParams1);
    [UincParams2,Th2] = Scat2TypeShifted.UincOnFieldHelper(ExtPrb2,ExParams2);
    
    u_tot1 = u1 + Uinc(UincParams1,Th1,IncAng,k);
    u_tot2 = u2Shifted + Uinc(UincParams2,Th2,IncAng,k);
    %u_tot3 = u3 + Uinc(UincParams2,Th2,IncAng,k);
    u_tot3 = u3 + Uinc(UincParams1,Th1,IncAng,k);
    
    
    fig=figure;
    subplot(131), DrawField(u_tot1,ExtPrb1,Setup1,ExParams1,IncAng,'Tot field1',fig);
    subplot(132), DrawField(u_tot2,ExtPrb2,Setup2Shited,ExParams2,IncAng,'Tot field2(shifted)',fig);
    subplot(133), DrawField(u_tot3,ExtPrb1,Setup1,ExParams1,IncAng,'Tot field3(f1 shifted)',fig);
    
    %figure
    %Quiver(u1,u3,ExtPrb1,Setup1,ExParams1,IncAng,'Tot field1',fig);
end

[FFP_plus2,FFP_minus2] = getFFP(u2Shifted,Setup2Shited,ExParams2);
[FFP_plus3,FFP_minus3] = getFFP(u3,Setup1,ExParams1);

%%
figure('units', 'normalized', 'position', [0.1 0.1 0.8, 0.4])

phi = FFP_minus3.phi_refl;
xhat = [cos(phi),sin(phi)];
% IncAng = IncAng-pi;
KincHat = [cos(IncAng),sin(IncAng)];
arg = (xhat(1)-KincHat(1))*Shift(1) + (xhat(2)-KincHat(2))*Shift(2); Title = 'FFP minus from 10th last circle';
%     arg =- arg; title('FFP plus, -arg'), Title = 'FFP plus, -arg'; warning('wrong sign!')


ShitftedF0 = FFP_minus3.ffp; %approx using
F2 = FFP_minus2.ffp; % computed from the field

sgtitle(Title);
subplot(221)
plot(phi*180/pi,real(F2) ,'bo',phi*180/pi,real(ShitftedF0),'r+'), title('real')
subplot(222)
plot(phi*180/pi,imag(F2) ,'bo',phi*180/pi,imag(ShitftedF0),'r+'), title('imag')
subplot(223)
plot(phi*180/pi,abs(F2) ,'bo',phi*180/pi,abs(ShitftedF0),'r+'), title('abs')
subplot(224)
plot(phi*180/pi,angle(F2),'bo',phi*180/pi,angle(ShitftedF0),'r+'), title('angle')

diff = F2 - ShitftedF0;


figure('units', 'normalized', 'position', [0.1 0.1 0.8, 0.4])
sgtitle(Title);
subplot(131)
plot(phi*180/pi,real(diff)), title('real')
subplot(132)
plot(phi*180/pi,imag(diff)), title('imag')
subplot(133)
plot(phi*180/pi,angle(diff)), title('angle')


fprintf('FFP minus: max-err(F1-ShiftedF0)=%f, max-err(angle(F1-ShiftedF0))=%f\n', norm(diff,inf)/norm(F2,inf), norm(angle(diff),inf))

if 0
    %%
    figure('units', 'normalized', 'position', [0.1 0.1 0.8, 0.4])
    
    phi = FFP_plus1.phi;
    xhat = [cos(phi),sin(phi)];
    IncAng = IncAng+pi;
    KincHat = [cos(IncAng),sin(IncAng)];
    arg = (xhat(1)-KincHat(1))*Shift(1) + (xhat(2)-KincHat(2))*Shift(2); Title = 'FFP plus, arg';
    %argwrong =- arg; title('FFP plus, -arg'), Title = 'FFP plus, -arg';
    
    
    F0 = FFP_plus1.ffp;
    ShitftedF0 = F0.*exp(1i*k*arg);
    %ShitftedF0wrong = F0.*exp(1i*k*argwrong);
    F1 = FFP_plus2.ffp;
    
    %diff = F1 - ShitftedF0;
    %diffwrong = F1 - ShitftedF0wrong;
    %plot(phi*180/pi,real(diff),'b',phi*180/pi,-real(diffwrong),'r')
    
    sgtitle(Title);
    subplot(131)
    plot(phi*180/pi,real(F1),'bo',phi*180/pi,real(ShitftedF0),'r+'), title('real')
    subplot(132)
    plot(phi*180/pi,imag(F1),'bo',phi*180/pi,imag(ShitftedF0),'r+'), title('imag')
    subplot(133)
    plot(phi*180/pi,angle(F1),'bo',phi*180/pi,angle(ShitftedF0),'r+'), title('angle')
    %
    %
    %     fprintf('FFP plus: max-err(F1-ShiftedF0)=%f, max-err(angle(F1-ShiftedF0))=%f\n', norm(diff,inf), norm(angle(diff),inf))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,ExtPrb,cn0,cn1] = Solve(BC,Setup)
    
    ExtPrb =  Solvers.ExteriorSolver( Setup );
    
    xi = zeros(Setup.Grid.Nr,Setup.Grid.Nth);%spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
    [xi(ExtPrb.GridGamma),cn0,cn1] = ExtPrb.xi(BC,ExtPrb.GridGamma);
    
    if 0 % some old debugging test with known exact, not sure what it was
        bn=5; %should this be HandelIndex?
        xi(ExtPrb.GridGamma) = besselh(bn,HankelType,k*Grid.R(ExtPrb.GridGamma)).*exp(bn*1i*Grid.Theta(ExtPrb.GridGamma));
    end
    
    u = ExtPrb.P_Omega(xi);
end

function [FFP_plus,FFP_minus] = getFFP(u,Setup,ExParams)
    k = Setup.CoeffsParams.k;
    phi = linspace(0,2*pi,1000);
            
    if isa(u,'struct')
        field = u;
        theta = Setup.theta;
    else
        Grid = Setup.Grid;
        theta = Grid.theta;
        nn = numel(Grid.r) - 10;

        field.value = u(nn,:);
        field.normal_deriv = (8 * (u(nn+1,:) - u(nn-1,:)) - (u(nn+2,:) - u(nn-2,:)) ) / (12 * Grid.dx);
    end
    
    S = ScattererShape(ExParams,theta);
    curve.z = [S.x;S.y];
    curve.nz = [S.dx;S.dy];
    
    FFP_plus  = SAR.plus_i.SARUtils.doFFP_curve(k,phi,curve,field);
    FFP_minus = SAR.minus_i.mi__SARUtils.doFFP_curve(k,phi,curve,field);
end

function DrawFFP(u,Setup,ExParams)
    
    k = Setup.CoeffsParams.k;
    phi = linspace(0,2*pi,1000);
    
    Grid = Setup.Grid;
    nn = numel(Grid.r) - 10;
    
    
    S = ScattererShape(ExParams,Grid.theta);
    
    curve.z = [S.x;S.y];
    curve.nz = [S.dx;S.dy];
    
    
    field.value = u(nn,:);
    field.normal_deriv = (8 * (u(nn+1,:) - u(nn-1,:)) - (u(nn+2,:) - u(nn-2,:)) ) / (12 * Grid.dx);
    
    FFP_plus  = SAR.plus_i.SARUtils.doFFP_curve(k,phi,curve,field);
    FFP_minus = SAR.minus_i.mi__SARUtils.doFFP_curve(k,phi,curve,field);
    
    subplot(2,2,3)%[3,4])
    plot(FFP_plus.phi*180/pi,abs(FFP_plus.ffp))
    
    subplot(2,2,4)
    plot(FFP_minus.phi_refl*180/pi,abs(FFP_minus.ffp))
    
end

function DrawField(u,Prb,Setup,ExParams,IncAng, TitleMsg,fig)
    
    Grid = Setup.Grid;
    k = Setup.CoeffsParams.k;
    
    R  = ones(size(Grid.R)).*NaN;
    Th = ones(size(Grid.R)).*NaN;
    Nm = Prb.Scatterer.Nm;
    R(Nm) = Grid.R(Nm);
    Th(Nm)= Grid.Theta(Nm);
    %R = Grid.R;
    %Th= Grid.Theta;
    
    R(:,Grid.Nth+1)=R(:,Grid.Nth);
    Th(:,Grid.Nth+1)=2*pi;
    
    X = (R .* cos(Th));
    Y = (R .* sin(Th));
    
    U= ones(size(Grid.R)).*NaN;
    U(Nm)=u(Nm);
    %U=u;
    U(:,Grid.Nth+1)=u(:,1);
    
    
    mypcolor(X,Y,abs(U) , [TitleMsg ', abs, inc ang=' num2str(IncAng*180/pi) 'degree' ], fig, @() DrawScatterrerShape(ExParams,Setup));
    
    
    %     UincParams1 = ExParams;
    %     UincParams1.r = Prb.Scatterer.R;
    %     uinc = Uinc(UincParams1,Prb.Scatterer.Th,IncAng,k);
    %
    %
    %     subplot(2,2,2)
    %     mypcolor(X,Y,abs(uinc+U) , ['total field, abs, inc ang=' num2str(IncAng*180/pi) 'degree' ], fig, @() DrawScatterrerShape(ExParams));
end

function S = ScattererShape(Params,Angle)
    
    switch Params.ScattererType
        case Tools.Enums.Scatterer.Circle
            S = struct('x',Params.r*cos(Angle),'y',Params.r*sin(Angle), ...
                'dx',-Params.r*sin(Angle),'dy',-Params.r*cos(Angle));
            
        case Tools.Enums.Scatterer.Ellipse
            S = struct('x',Params.FocalDistance*cosh(Params.eta).*cos(Angle),'y', Params.FocalDistance*sinh(Params.eta).*sin(Angle), ...
                'dx',-Params.FocalDistance*cosh(Params.eta).*sin(Angle),'dy', Params.FocalDistance*sinh(Params.eta).*cos(Angle));
            
        case Tools.Enums.Scatterer.StarShaped
            [x,dx] = Params.Parameterization.XHandle.Derivatives(Angle);
            [y,dy] = Params.Parameterization.YHandle.Derivatives(Angle);
            S = struct('x',x,'y',y,'dx',dx,'dy',dy);
            
    end
    
end

function DrawScatterrerShape(Params,Setup)
    th=0:0.0001:2*pi;
    
    S = ScattererShape(Params,th);
    
    plot(S.x,S.y,'k','LineWidth',2)
    
    if nargin>1
        r0 = Setup.Grid.x(1);
        r1 = Setup.Grid.x(end);
        
        plot(r0*cos(th),r0*sin(th),'b--','LineWidth',2)
        plot(r1*cos(th),r1*sin(th),'b--','LineWidth',2)
    end
    %     switch Params.ScattererType
    %         case Tools.Enums.Scatterer.Circle
    %             plot(Params.r*cos(th),Params.r*sin(th),'k','LineWidth',2)
    %
    %         case Tools.Enums.Scatterer.Ellipse
    %             plot(Params.FocalDistance*cosh(Params.eta).*cos(th),Params.FocalDistance*sinh(Params.eta).*sin(th),'k','LineWidth',2)
    %
    %         case Tools.Enums.Scatterer.StarShaped
    %             plot(Params.Parameterization.XHandle.Derivatives(th),Params.Parameterization.YHandle.Derivatives(th),'k','LineWidth',2)
    %
    %     end
    
end


function uinc = Uinc(Params,phi,IncAng,k)
    
    switch Params.ScattererType
        case Tools.Enums.Scatterer.Circle
            x = Params.r .* cos(phi);
            y = Params.r .* sin(phi);
            
        case Tools.Enums.Scatterer.Ellipse
            x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
            y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
            
        case Tools.Enums.Scatterer.StarShaped
            
            try
                x = Params.Parameterization.XHandle.Derivatives(phi);
                y = Params.Parameterization.YHandle.Derivatives(phi);
            catch
                x = Params.r.*cos(phi);
                y = Params.r.*sin(phi);
            end
            
    end
    
    uinc = exp( 1i.* k .* (x.*cos(IncAng) + y.*sin(IncAng)) );
end

function duinc = detaUinc(Params,phi,IncAng,k)
    
    switch Params.ScattererType
        case Tools.Enums.Scatterer.Circle
            dx = cos(phi);
            dy = sin(phi);
            
        case Tools.Enums.Scatterer.Ellipse
            dx = Params.FocalDistance * sinh(Params.eta) .* cos(phi);
            dy = Params.FocalDistance * cosh(Params.eta) .* sin(phi);
            
        case Tools.Enums.Scatterer.StarShaped
            [x,dx] = Params.Parameterization.XHandle.Derivatives(phi);
            [y,dy] = Params.Parameterization.YHandle.Derivatives(phi);
            
    end
    
    uinc = Uinc(Params,phi,IncAng,k)  ;
    
    duinc = 1i .* k .*  uinc .* (dx.*cos(IncAng) + dy.*sin(IncAng));
    %   h = FocalDist*sqrt(sinh(eta).^2 + sin(phi).^2);
    % duinc = duinc./h;
    
    if Params.ScattererType == Tools.Enums.Scatterer.StarShaped
        h = sqrt(dx.^2 + dy.^2);
        duinc = 1i .* k .*  uinc .* (dy.*cos(IncAng) - dx.*sin(IncAng))./h;
    end
    
end


function Quiver(u,v,Prb,Setup,ExParams,IncAng, TitleMsg,fig)
    
    Grid = Setup.Grid;
    k = Setup.CoeffsParams.k;
    
    R  = ones(size(Grid.R)).*NaN;
    Th = ones(size(Grid.R)).*NaN;
    Nm = Prb.Scatterer.Nm;
    R(Nm) = Grid.R(Nm);
    Th(Nm)= Grid.Theta(Nm);
    %R = Grid.R;
    %Th= Grid.Theta;
    
    R(:,Grid.Nth+1)=R(:,Grid.Nth);
    Th(:,Grid.Nth+1)=2*pi;
    
    X = (R .* cos(Th));
    Y = (R .* sin(Th));
    
    U= ones(size(Grid.R)).*NaN;
    U(Nm)=u(Nm);
    %U=u;
    U(:,Grid.Nth+1)=u(:,1);
    
    V= ones(size(Grid.R)).*NaN;
    V(Nm)=v(Nm);
    %U=u;
    V(:,Grid.Nth+1)=v(:,1);
    
    
   % mypcolor(X,Y,abs(U) , [TitleMsg ', abs, inc ang=' num2str(IncAng*180/pi) 'degree' ], fig, @() DrawScatterrerShape(ExParams,Setup));
    quiver(X,Y,U,V)
    
    %     UincParams1 = ExParams;
    %     UincParams1.r = Prb.Scatterer.R;
    %     uinc = Uinc(UincParams1,Prb.Scatterer.Th,IncAng,k);
    %
    %
    %     subplot(2,2,2)
    %     mypcolor(X,Y,abs(uinc+U) , ['total field, abs, inc ang=' num2str(IncAng*180/pi) 'degree' ], fig, @() DrawScatterrerShape(ExParams));
end

function u = Assemble(theta, BasisHandle,Indices,coeffs)
    BasisAddParams=[]; %required for Chebyshev Basis, that is not supposed to be used here
    u = 0;
    for j=1:numel(Indices)
        uj = BasisHandle(theta,Indices(j),BasisAddParams);
        u = u + coeffs(j).*uj.xi0;
    end    
end
