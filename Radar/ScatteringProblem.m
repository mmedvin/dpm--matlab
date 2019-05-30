function u = ScatteringProblem(WaveNumber,IncidentAngle,BoundaryConditions,ScatterrerType)
    
    if(nargin ==0)
        [k,IncAng,BC,ScatType] = DefaultInput();
    else
        [k,IncAng,BC,ScatType] = deal(WaveNumber,IncidentAngle,BoundaryConditions,ScatterrerType);
    end
    
    LastKnownVersionOfTools = '2.0.0.0';
    Tools.VerifyVersion(LastKnownVersionOfTools);
    Ngrids = 3;
    
    %Ellipse Axes:
    a=1;
    b=a/2;
    
    %Ring radiuses:
    r0=0.3; r1=2.2; %r0=0.8; r1=2.2;
    %r0=0.8; r1=2.2;
    
    KindOfConvergance = Tools.Enums.Convergance.Grid; %Exact; %
    
    BasisType = Tools.Enums.Basis.Fourier;
    ChebyshevRange = struct('a',-pi,'b',pi);%parameter for Chebyshev  Basis if used, don't change it!
    
    %special test with known exact solution
    HankOrPlane = 'PlaneWave';% 'PlaneWave' or 'Hankel'
    HankelIndex = 3;
    HankelType = 2;
    
    %Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b)); %needed for test only
    %Parameterization  = Tools.Parameterizations.ParametricKite(struct('a',1,'b',.65*2,'c',1.5));
    Parameterization  = Tools.Parameterizations.ParametricSubmarine(struct('a',1.8,'b',1.8/3,'c',2,'p',100));
    
    %shift=[1/3,1/3];        
   %Parameterization  = Tools.Parameterizations.ParametricEllipse2(struct('a',a,'b',b,'xcenter',0,'ycenter',0.1,'rotation',0));

    
    
    [ScattererHandle,ScattererParams,Extension,ExParams] = ScatType.Helper(a,b,HankelIndex,HankelType,Parameterization);
    
    fprintf('Scattering Problem about %s\n',ScatType);
    ScatType.Print(a,b,HankelIndex,HankelType,Parameterization);
    fprintf('Polar Grid in Ring between r0=%-5.2f and r1=%-5.2f \n ',r0,r1 );
    
    
    if strcmpi(KindOfConvergance,'Exact') && strcmpi(HankOrPlane,'PlaneWave') && ~(a==1 && b==0.8 && r0 == 0.7*b && r1 == 1.8*a)
        error('are you sure?')
    end
    
    fprintf('BC:%s, Convergance:%s, data is %s, scatterer is %s, Basis is %s, Hnkl_n=%d \n',...
        BC.toString,KindOfConvergance.toString,HankOrPlane,ScatType,BasisType.toString,HankelIndex);
    
    
    if strcmpi(HankOrPlane,'PlaneWave')
        f1      = @(phi) -Uinc(ExParams,phi,IncAng,k);
        dfdn    = @(phi) -detaUinc(ExParams,phi,IncAng,k);
    elseif strcmpi(HankOrPlane,'Hankel')
        f1      = @(phi) ExactHank(ExParams,phi,k);
        dfdn    = @(phi) dnExactHank(ExParams,phi,k);
    end
    
    Basis = BasisType.Helper(f1,dfdn,[1e-06,1e-06]);
    
    Setup = struct(...
        'Basis'           , Basis, ...
        'Grid'            , Tools.Grid.PolarGrids(r0,r1,512,512), ...
        'CoeffsHandle'    , @Tools.Coeffs.ConstantWaveNumber, ...
        'CoeffsParams'    , struct('k',k,'r0',1.6), ...
        'ScattererHandle' , ScattererHandle, ...
        'ScattererParams' , ScattererParams, ...
        'CollectRhs'      , 0, ... %i.e. not
        'Extension'       , Extension, ...
        'ExtensionParams' , [] ...
        );
    
    [u,Prb,Setup] = ConverganceTest(Ngrids,BC,KindOfConvergance,Setup);
    
    fig = figure;
    
    subplot(2,2,1)
    DrawField(u,Prb,Setup,ExParams,IncAng,'scatterred field',fig);

    subplot(2,2,2)
    [UincParams,Th] = ScatType.UincOnFieldHelper(Prb,ExParams);
    u_tot = u + Uinc(UincParams,Th,IncAng,k);
    DrawField(u_tot,Prb,Setup,ExParams,IncAng,'total field',fig);

    
     DrawFFP(u,Setup,ExParams)
    
     saveas(gcf,[ScatType.toString(), num2str(IncAng*180/pi) '.jpg'],'jpg')
     
    if nargout==0, clear u, end
    
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


function [u,ExtPrb,Setup] = ConverganceTest(Nmax,BC,KindOfConvergance,Setup)
    
    k = Setup.CoeffsParams.k;
    r0 = Setup.Grid.x(1);
    r1 = Setup.Grid.x(end);
    
    ErrPre = 0;
    etinf=[];
    
    if KindOfConvergance == Tools.Enums.Convergance.Grid, Nmax = Nmax+1; end
    
    for n=1:Nmax %run different grids
        tic
        p=6;%1;
        Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
        
        Setup.Grid = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
        [u,ExtPrb] = Solve(BC,Setup);
        
        t=toc;
        
        % % % % % % % % % % % % % % % %
        % Comparison
        % % % % % % % % % % % % % % % %
        if strcmpi(KindOfConvergance,'Grid')
            if n > 1
                u1= spalloc(Nr,Nth-1,nnz(u0));
                u1(ExtPrb.Nm) = u0(ExtPrb.Nm);
                
                tmp = u(1:2:end,1:2:end)-u1(1:2:end,1:2:end);
                
                %etinf(n) =norm(tmp(:),inf);
                ErrTot = norm(tmp(:),inf);
                %fprintf('k=%d,M=%d,N=%-10dx%-10d\t etinf=%d\ttime=%d\n',k,Basis.M, Nr,Nth,full(etinf(n)),t);
                fprintf('k=%d,NBss0=%d,NBss1=%d,N=%-10dx%-10d\t ErrTot=%d\t rate=%-5.2f\t time=%d\n',k,Setup.Basis.NBss0,Setup.Basis.NBss1, Nr,Nth,ErrTot,log2(ErrPre/ErrTot),t);
                
                ErrPre = ErrTot;
                
                
                % fprintf('k=%-5.4f,pf=%-5.4f,pf=%-5.4f,IncAng=%d,M=%d,N=%-10dx%-10d\t etinf=%d\ttime=%d\t \n',k,  (k^5)* Grid.dx^4, (k^5) * Grid.dy^4,IncAngD,Basis.M, Nr,Nth,full(etinf(n)),t);
                
            end
            
            u0=spalloc(Nr*2-1,Nth*2-2,nnz(u));
            u0(1:2:end,1:2:end)=u;
            
        elseif strcmpi(KindOfConvergance,'Exact')
            %            phi = linspace(0,2*pi,300);
            phi = linspace(-pi,pi,300);
            
            app = 0;
            
            if strcmpi(Problem , 'Dirichlet')
                for j=1:numel(Basis.Indices)
                    bj = Basis.Handle(phi,Basis.Indices(j),Basis.AddParams);
                    app = app + cn1(j).*bj.xi0;
                end
                
                if strcmpi(HankOrPlane,'PlaneWave')
                    ex = ell_du_dn_exact_SM(a,b,k,phi,IncAng)';
                elseif strcmpi(HankOrPlane,'Hankel')
                    %                         ex = dnExactHank(ExParams,phi,k);????
                    ex = dnExactHank(ExParams,-phi,k);
                end
                
                % error('Dirichlet problem not implemented yet')
            elseif strcmpi(Problem , 'Neumann')
                for j=1:numel(Basis.Indices)
                    bj = Basis.Handle(phi,Basis.Indices(j), ExtPrb.Scatterer.MetricsAtScatterer,Basis.AddParams);
                    app = app + cn0(j).*bj.xi0;
                    %                 for j=-M:M
                    %                     app= app +  cn0(j+M+1).*exp(1i*j*phi);
                    %du_app = du_app + cn1(j+M+1).*exp(1i*j*phi);
                end
                
                if strcmpi(HankOrPlane,'PlaneWave')
                    ex = ell_hard_exact_Sol(a,b,k,IncAng,phi)';
                elseif strcmpi(HankOrPlane,'Hankel')
                    ex = ExactHank(ExParams,phi,k);
                end
            else
                error('Solving only Dirichlet or Neumann problem')
            end
            
            ex=ex/norm(ex,inf);
            app=app/norm(app,inf);
            
            tmp = app(:) - ex(:);
            etinf(n) =norm(tmp(:),inf);
            fprintf('k=%d,IncAng=%d,M=%d,N=%-10dx%-10d\t etinf=%d\ttime=%d\n',k,IncAngD,Basis.M, Nr,Nth,full(etinf(n)),t);
            
            
        else
            error('undefinded comparison method');
        end
    end
    
    fprintf('\n');
    
    if nargout==0
        clear u
    end
end

function [u,ExtPrb] = Solve(BC,Setup)
    
    ExtPrb =  Solvers.ExteriorSolver( Setup );
    
    xi = zeros(Setup.Grid.Nr,Setup.Grid.Nth);%spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
    xi(ExtPrb.GridGamma) = ExtPrb.xi(BC,ExtPrb.GridGamma);
    
    if 0 % some old debugging test with known exact, not sure what it was
        bn=5; %should this be HandelIndex?
        xi(ExtPrb.GridGamma) = besselh(bn,HankelType,k*Grid.R(ExtPrb.GridGamma)).*exp(bn*1i*Grid.Theta(ExtPrb.GridGamma));
    end
    
    u = ExtPrb.P_Omega(xi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function u = ExactHank(Params,phi,k)
    if strcmpi(Params.ScattererType,'ellipse')
        x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
        y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        x = Params.r .* cos(phi);
        y = Params.r .* sin(phi);
        %         r= Params.r;
        %         th = phi;
    end
    
    Z=x+1i*y;
    r=abs(Z);
    th=angle(Z);
    
    
    
    bn=Params.HankelIndex;
    u = besselh(bn,Params.HankelType,k*r).*exp(bn*1i*th);
end

function un = dnExactHank(Params,phi,k)
    
    if strcmpi(Params.ScattererType,'ellipse')
        x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
        y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        x = Params.r .* cos(phi);
        y = Params.r .* sin(phi);
        %         r= Params.r;
        %         th = phi;
    end
    
    Z=x+1i*y;
    r=abs(Z);
    th=angle(Z);
    
    bn=Params.HankelIndex;
    
    dudr  = k*0.5*(besselh(bn-1,Params.HankelType,k*r) - besselh(bn+1,Params.HankelType,k*r)).*exp(bn*1i*th);
    
    if strcmpi(Params.ScattererType,'ellipse')
        
        dudth  = 1i*bn*besselh(bn,Params.HankelType,k*r).*exp(bn*1i*th);
        
        drdeta = Params.FocalDistance^2 * cosh(Params.eta).*sinh(Params.eta) ./r;
        dthdeta = Params.FocalDistance^2 * cos(phi).*sin(phi) ./(r.^2);
        un = dudr.*drdeta + dudth.*dthdeta;
        
        %         SM = EllipticalMetrics(Params.FocalDistance,Params.eta,phi);
        % un = un.*SM.h;
    elseif strcmpi(Params.ScattererType,'circle')
        un = dudr;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function [k,IncAng,BC,ScatType] = DefaultInput()
    k           = 5;
    IncAng      = pi/5;
    BC          = Tools.Enums.BoundaryConditions.Dirichlet;%;Neumann;%
    ScatType    = Tools.Enums.Scatterer.StarShaped;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
