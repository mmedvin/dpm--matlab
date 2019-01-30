function SARsimulation(filename,AR,GridParam,ScattererType, K0, Phi)
% clear all, close all, clc
    if nargin == 0
        GridParam=7;
        dPhi=0.04;%0.004;
        dK = 1;%0.5;
        Phi = pi*(-0.2:dPhi:0.2);%0;%-10:5:10;% 
        K0 = 50:dK:55;%5;%5:7;%50;%        
        ScattererType = 'ellipse'; %'circle';%'star' , 'ellipse'
        AR=2;
        filename = [pwd filesep 'SAR.mat'];
    end
    
   
    
       Solve4AllKInc(K0,Phi,AR,GridParam,ScattererType,@Uinc, filename);
%     PlotM(filename);
    
%     x=linspace(-1,1,1000);
%     y=x;
%      GenerateSARImage(x,y,filename)
 %   load(filename);

 %test
%  th = linspace(0,2*pi,100)
%  uinf = zeros(size(th));
%  IncAng = 33 * pi/180;
%  for i=1:numel(th)
%      s=th(i);
%      uinf(i) = FFP(1,1,5, IncAng,cn, Parameterization, Basis)
%  end
end   
    
function Solve4AllKInc(K0,Ang0,AR,GridParam,ScattererType,Uinc,filename)
     t1=tic;

    NHR = 1.6; %variable wavenumber param - do not change
    
   %size of problem
    Nr=2^GridParam;	Nth=2^GridParam;
    Nx=2^GridParam;	Ny=2^GridParam;

%     Z = zeros(Nr,Nth);
%     SD=struct('Extu',Z,'cn0',[],'cn1',[],'k',0,'phi',0);
%     M=struct('scattData',  repmat({SD},numel(K0),numel(Ang0)));
    
    Basis   = Tools.Basis.FourierBasis.BasisHelper(@sin,@sin,[100,100]);
    
    WaveNumberHandle = @Tools.Coeffs.ConstantWaveNumber;
 
    if  strcmpi(ScattererType,'ellipse')
        ScattererType = 'StarShapedScatterer';
        
        a=1;%b=a/2;
        b=a/AR;
        %exterior problem in ring
        r0 = 0.7*b;
        r1 = 1.8*a; %ABC set on that circle
        
        %interior problem in square
        x1=-(a+0.2);xn=(a+0.2);
        y1=-(b+0.2);yn=(b+0.2);
        
        FocalDistance = sqrt(a^2-b^2);
        Eta0 = acosh(a/FocalDistance);
        
        Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b));
        
        ScattererHandle  = @Tools.Scatterer.StarShapedScatterer;
        UParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization);

        ScattererParams  = UParams;
        ScattererParams.Stencil=9;
        
        Extension =  @Tools.Extensions.TwoTupleExtension;

    elseif  strcmpi(ScattererType,'star')
        
        ScattererType = 'StarShapedScatterer';
        ScattererHandle  = @Tools.Scatterer.StarShapedScatterer;
        Parameterization  = Tools.Parameterizations.ParametricStar();
                
        x1=-1.7;xn=1.7;
        y1=-1.7;yn=1.7;
        
        r0 = 0.8;
        r1 = 2.2; %ABC set on that circle

        UParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization);

        ScattererParams  = UParams;
        ScattererParams.Stencil=9;
        
        Extension =  @Tools.Extensions.TwoTupleExtension;
    
    elseif strcmpi(ScattererType,'circle')
        ScattererHandle  = @Tools.Scatterer.PolarScatterer;
        
        R0 =1; %radius of scatterrer
        
        x1=-1.2;xn=1.2;
        y1=-1.2;yn=1.2;
        
        r0 = 0.8;
        r1 = 2.0; %ABC set on that circle
      
        ScattererParams  = struct('r0',R0,'ExpansionType',15, 'Stencil', 9);
        UParams = struct('ScattererType','circle', 'r', R0);
        
        Extension =  @Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension;
    else
        error('wrong scatterrer')
    end
    
    PlrGrid = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
    CrtsGrid = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);

    % PsiI  = cell(numel(K0),numel(Ang0));
    % PsiII = cell(numel(K0),numel(Ang0));
    % uinc  = zeros(numel(K0),numel(Ang0));
    
    fprintf('\n SARsimulation, Ellipse a=%4.2f, b=%4.2f , filename=%s\n',a,b,filename);

    for indx=1:numel(K0) 
        k0=K0(indx);
        k1=1.1*k0;  
        
       
       % fprintf('starting k0=%5.3f,k1=%5.3f after %6.4f secs from the begining\n',k0,k1,toc(t1));
        
%  k1=10*k0;

        ExtPrb = Solvers.ExteriorSolver( struct('Basis',Basis, ...
            'Grid', PlrGrid, 'CoeffsHandle', WaveNumberHandle, ...
            'CoeffsParams', struct('k',k0,'r0',NHR), ...
            'ScattererHandle',ScattererHandle, ...
            'ScattererParams', ScattererParams, 'CollectRhs',0, ... %i.e. no
            'Extension', Extension, 'ExtensionParams',[] ));
                
        IntPrb = Solvers.InteriorHomoSolver( struct('Basis',Basis, ...
            'Grid', CrtsGrid, 'CoeffsHandle', WaveNumberHandle, ...
            'CoeffsParams', struct('k',k1,'r0',NHR), ...
            'ScattererHandle',ScattererHandle, ...
            'ScattererParams', ScattererParams, ...
            'CollectRhs',1, ... %i.e. yes
            'Extension', Extension, 'ExtensionParams',[] ));

        %t2=toc(t1);
        %fprintf('computing diff inc angles for k0=%5.3f,k1=%5.3f after %6.4f secs\n',k0,k1,t2);
        
        for jndx = 1:numel(Ang0)            
                        
            IncAng = Ang0(jndx);%*pi/180;
            
            UincParams = UParams;
            UincParams.r = ExtPrb.Scatterer.r;
            %UincParams  = UParams;struct('ScattererType','circle','r',ExtPrb.Scatterer.r);
            
            rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
            uinc = Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k0);
            rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
            
            cn = [ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ] \ rhs;
   
            if 1
                Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = [IntPrb.W{1}(IntPrb.GridGamma,:),IntPrb.W{2}(IntPrb.GridGamma,:)]*cn;
                Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}(ExtPrb.GridGamma,:),ExtPrb.W{2}(ExtPrb.GridGamma,:)]*cn ;%- uinc;
                
                UincParams = UParams;
                UincParams.r = PlrGrid.R;
                %UincParams  = struct('ScattererType','circle','r',PlrGrid.R);
                uinc = Uinc(UincParams,PlrGrid.Theta,IncAng,k0);
                
                Extu = ExtPrb.P_Omega(Extxi,uinc );
                
                %create SAR image
                if  strcmpi(ScattererType,'circle')
                    s= ExtPrb.Scatterer.th;
                else
                    s= ExtPrb.Scatterer.nrml_t;
                end
                %            uinf(indx,jndx) = FFP(alpha,s,k0, IncAng, cn,Parameterization, Basis);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             R = PlrGrid.R;
                %             Th = PlrGrid.Theta;
                %
                %             R(:,Nth)=R(:,Nth-1);
                %             Th(:,Nth)=2*pi;
                %             nExtu = Extu;
                %             nExtu(:,Nth)=nExtu(:,1);
                %
                %             All=Intu;
                %             All(IntPrb.Scatterer.Mm) = griddata(R .* cos(Th),R .* sin(Th) ,full(nExtu), CrtsGrid.X(IntPrb.Scatterer.Mm), CrtsGrid.Y(IntPrb.Scatterer.Mm) );
                %
                %             D{indx,jndx} = All;
                M.scattData{indx,jndx}.Extu = Extu;
                %M{indx,jndx}.Intu = Intu;
            end
            
            %save it all
            M.scattData{indx,jndx}.cn0 = cn(1:Basis.NBss0);
            M.scattData{indx,jndx}.cn1 = cn( (Basis.NBss0+1):end);
            M.scattData{indx,jndx}.k=k0;
            M.scattData{indx,jndx}.phi = IncAng;
        end
               
        t3=toc(t1);
        fprintf('done computing diff inc angles for k0=%5.3f,k1=%5.3f after %6.4f secs\n',k0,k1,toc(t1));
    end

    
    M.Basis = Basis;
    M.PlrGrid=PlrGrid;

    save(filename);%, 'PlrGrid', 'CrtsGrid' ,'ExtPrb', 'IntPrb');

    fprintf(' saving took %6.4f secs\n',toc(t1)-t3);

    
     fprintf('finished all after %6.4f secs in total\n',toc(t1));

    
    %  save('SAR','D','SarIm', 'PsiI', 'PsiII','uinf','CrtsGrid','PlrGrid');

%     save('SAR','M','SarIm', 'uinf','CrtsGrid','PlrGrid');
%     mypcolor(x,x,abs(SarIm));
        
        %%%%%%%%%%%%%%%
          
end


function [uinc,uinc_t,uinc_tt,uinc_3t,uinc_4t] = Uinc(Params,phi,IncAng,k)
    
    IsStarshaped = false;
    
    if strcmpi(Params.ScattererType,'ellipse')
        x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
        y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        x = Params.r .* cos(phi);
        y = Params.r .* sin(phi);
    elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
        IsStarshaped = true;
        try
            x = Params.Parameterization.XHandle.Derivatives(phi);
            y = Params.Parameterization.YHandle.Derivatives(phi);
        catch
            x = Params.r.*cos(phi);
            y = Params.r.*sin(phi);
        end        
    end
    
    uinc = exp( 1i.* k .* (x.*cos(IncAng) + y.*sin(IncAng)) );
    
    if nargout > 1 && IsStarshaped
        %try
        [x,xt,xtt,x3t,x4t] = Params.Parameterization.XHandle.Derivatives(phi);
        [y,yt,ytt,y3t,y4t] = Params.Parameterization.YHandle.Derivatives(phi);
        %catch
        %x = Params.r.*cos(phi);
        %y = Params.r.*sin(phi);
        %end
        
        uinc_t  = 1i .* k .*  uinc    .* (xt .*cos(IncAng) + yt .*sin(IncAng));
        uinc_tt = 1i .* k .* (uinc_t  .* (xt .*cos(IncAng) + yt .*sin(IncAng)) +     uinc    .* (xtt.*cos(IncAng) + ytt.*sin(IncAng)) );
        uinc_3t = 1i .* k .* (uinc_tt .* (xt .*cos(IncAng) + yt .*sin(IncAng)) + 2 * uinc_t  .* (xtt.*cos(IncAng) + ytt.*sin(IncAng)) + uinc .* (x3t.*cos(IncAng) + y3t.*sin(IncAng)));
        uinc_4t = 1i .* k .* (uinc_3t .* (xt .*cos(IncAng) + yt .*sin(IncAng)) + 3 * uinc_tt .* (xtt.*cos(IncAng) + ytt.*sin(IncAng)) ...
            +       3  *  uinc_t  .* (x3t.*cos(IncAng) + y3t.*sin(IncAng)) +     uinc    .* (x4t.*cos(IncAng) + y4t.*sin(IncAng)));
    end
    
end


function [duinc,duinc_t,duinc_tt] = detaUinc(Params,phi,IncAng,k)
    IsStarshaped = false;
    if strcmpi(Params.ScattererType,'ellipse')
        dx = Params.FocalDistance * sinh(Params.eta) .* cos(phi);
        dy = Params.FocalDistance * cosh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        dx = cos(phi);
        dy = sin(phi);
    elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
        IsStarshaped = true;
        [x,dx] = Params.Parameterization.XHandle.Derivatives(phi);
        [y,dy] = Params.Parameterization.YHandle.Derivatives(phi);
    end
    
    uinc = Uinc(Params,phi,IncAng,k)  ;
    
    duinc = 1i .* k .*  uinc .* (dx.*cos(IncAng) + dy.*sin(IncAng));
    %   h = FocalDist*sqrt(sinh(eta).^2 + sin(phi).^2);
    % duinc = duinc./h;
    
    if IsStarshaped
        
        h = sqrt(dx.^2 + dy.^2);
        duinc = 1i .* k .*  uinc .* (dy.*cos(IncAng) - dx.*sin(IncAng))./h;
        
        if nargout > 1
            [x,xt,xtt,x3t,x4t] = Params.Parameterization.XHandle.Derivatives(phi);
            [y,yt,ytt,y3t,y4t] = Params.Parameterization.YHandle.Derivatives(phi);
            
            [uinc,uinc_t,uinc_tt] = Uinc(Params,phi,IncAng,k);
            
            ht  = (xt.*xtt + yt.*ytt)./h;
            htt = (xtt.^2 + ytt.^2 + xt.*x3t + yt.*y3t - ht.^2)./h;
            h3t = (3*xtt.*x3t + 3*ytt.*y3t + xt.*x4t + yt.*y4t +  - 3*ht.*htt)./h;
            
            duinc_t = 1i .* k .* ( uinc_t  .* (yt .*cos(IncAng) - xt .*sin(IncAng))./h  + uinc   .* (ytt.*cos(IncAng) - xtt.*sin(IncAng))./h + uinc   .* (yt .*cos(IncAng) - xt .*sin(IncAng)).*(-ht./(h.^2)) );
            duinc_tt = 1i .* k .*( uinc_tt .* (yt .*cos(IncAng) - xt .*sin(IncAng))./h  + uinc_t .* (ytt.*cos(IncAng) - xtt.*sin(IncAng))./h + uinc_t .* (yt .*cos(IncAng) - xt .*sin(IncAng)).*(-ht./(h.^2)) ...
                +uinc_t  .* (ytt.*cos(IncAng) - xtt.*sin(IncAng))./h  + uinc   .* (y3t.*cos(IncAng) - x3t.*sin(IncAng))./h + uinc   .* (ytt.*cos(IncAng) - xtt.*sin(IncAng)).*(-ht./(h.^2)) ...
                +uinc_t  .* (yt .*cos(IncAng) - xt .*sin(IncAng)).*(-ht./(h.^2))  + uinc .* (ytt.*cos(IncAng) - xtt.*sin(IncAng)).*(-ht./(h.^2))  + uinc .* (yt.*cos(IncAng) - xt.*sin(IncAng)).*((2*ht.^2)./(h.^3) - htt./(h.^2)) ...
                );
            
            
            
            
        end
        
    end
    
end

function GenerateSARImage(x,y,filename)
    load(filename);

    alpha = exp(1i*pi/4)/sqrt(8*pi*k0);
    
    SarIm=zeros(numel(x), numel(y));
    for i=1:numel(x)
        for j=1:numel(x)
            for indx=1:numel(K0)
                k0=K0(indx);
                for jndx = 1:numel(Ang0)
                    IncAng = Ang0(jndx)*pi/180;
                    xhat=pi + IncAng;
                    SarIm(i,j)= uinf(indx,jndx)* exp(2*1i*k0.* (x(i).*cos(xhat) + y(j).*sin(xhat)) );
                    
                end
            end
        end
    end
    
    mypcolor(x,y,abs(SarIm));
end

function uinf = FFP(alpha,s,k0, IncAng,cn, Parameterization, Basis)
    xhat=pi + IncAng;
        
    if isempty(Parameterization)
    else
        [Zx,dZx] = Parameterization.XHandle.Derivatives(s);
        [Zy,dZy] = Parameterization.YHandle.Derivatives(s);
    end
    xz = (cos(xhat).*Zx + sin(xhat).*Zy);
    G = exp(-1i*k0.*xz);
    Gn = -1i .* k0 .*  G .* (dZx.*cos(xhat) + dZy.*sin(xhat));
    
    bn = Basis.Handle(s,Basis.Indices0);
    bn = bn.xi0.';
    PsiI = bn*Gn;
    
    bn = Basis.Handle(s,Basis.Indices1);
    bn = bn.xi0.';
    PsiII= bn*G;
    
    uinf = alpha* [PsiI.', -PsiII.']*cn;
end


function PlotM(filename)
    
    load(filename);%, 'PlrGrid', 'CrtsGrid' ,'ExtPrb', 'IntPrb');
    
    for indx=1:numel(K0)
        for jndx=1:numel(Ang0)
            
            %                     mypcolor(CrtsGrid.X, CrtsGrid.Y,abs(D{indx,jndx}),'total field, abs',[], @() draw_circle(R0,UParams));
            
            %                     %%%%%%%%%%%%%%%
            
            %                     R = PlrGrid.R;
            %                     Th = PlrGrid.Theta;
            %
            %                     R(:,Nth)=R(:,Nth-1);
            %                     Th(:,Nth)=2*pi;
            %                     nExtu = Extu;
            %                     nExtu(:,Nth)=nExtu(:,1);
            %
            %                     All=Intu;
            %                     All(IntPrb.Scatterer.Mm) = griddata(R .* cos(Th),R .* sin(Th) ,full(nExtu), CrtsGrid.X(IntPrb.Scatterer.Mm), CrtsGrid.Y(IntPrb.Scatterer.Mm) );
            %
            %                     mypcolor(CrtsGrid.X, CrtsGrid.Y,abs(All),'total field, abs',4, @() draw_circle(R0,UParams));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            R  = ones(size(PlrGrid.R)).*NaN;
            Th = ones(size(PlrGrid.R)).*NaN;
            Nm = ExtPrb.Scatterer.Nm;
            R(Nm) = PlrGrid.R(Nm);
            Th(Nm)= PlrGrid.Theta(Nm);
            
            R(:,Nth)=R(:,Nth-1);
            Th(:,Nth)=2*pi;
            
            XExt = (R .* cos(Th));
            YExt = (R .* sin(Th));
            
            tExtu= ones(size(PlrGrid.R)).*NaN;
            tExtu(Nm)=M{indx,jndx}.Extu(Nm);
            tExtu(:,Nth)=tExtu(:,1);
            
            Np=IntPrb.Scatterer.Np;
            tIntu = ones(size(CrtsGrid.X)).*NaN;
            tIntu(Np) = M{indx,jndx}.Intu(Np);
            %tIntu(IntPrb.Scatterer.Mm)=0;%NaN;
            
            XInt = ones(size(CrtsGrid.X)).*NaN;
            YInt = ones(size(CrtsGrid.X)).*NaN;
            
            XInt(Np) = CrtsGrid.X(Np);
            YInt(Np) = CrtsGrid.Y(Np);
            
            mypcolor([XInt,XExt],[YInt,YExt],abs([tIntu,tExtu]) , 'total field, abs' , [], @() draw_circle(UParams));
            %mypcolor([XInt,XExt],[YInt,YExt],real([tIntu,tExtu]), 'total field, real', 2, @() draw_circle(UParams));
            %mypcolor([XInt,XExt],[YInt,YExt],imag([tIntu,tExtu]), 'total field, imag', 3, @() draw_circle(UParams));
            
        end
    end   
end

function draw_circle(UParams)
    th=0:0.0001:2*pi;
    
    if  strcmpi(UParams.ScattererType,'StarShapedScatterer')
        plot(UParams.Parameterization.XHandle.Derivatives(th),UParams.Parameterization.YHandle.Derivatives(th),'k','LineWidth',2)
    elseif strcmpi(UParams.ScattererType,'circle')
        plot(UParams.r*cos(th),UParams.r*sin(th),'k','LineWidth',2)
    end
end



