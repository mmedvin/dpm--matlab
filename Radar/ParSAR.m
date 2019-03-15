function ParSAR(path,AR,GridParam,ScattererType, K0, Phi,rchoice,kfactor)
    % clear all, close all, clc
    if nargin == 0
        GridParam=7;
        dPhi=0.04;%0.004;
        dK = 2;%0.5;
        Phi = pi*(-0.2:dPhi:0.2);%0;%-10:5:10;%
        K0 = 50:dK:55;%5;%5:7;%50;%
        ScattererType = 'ellipse'; %'circle';%'star' , 'ellipse'
        AR=2;
        filename = [pwd filesep 'ParSAR.mat'];
        
        rchoice=[10,11,12];
        
        %spec.theta = pi * (0 : 0.0025 : 1.999);
        
        kfactor = 1.1;
        
    end
    
    filename = sprintf('ParSar_%sgp%dkf%g_',firstUpper(ScattererType), GridParam,kfactor);
    %filename = [path filesep filename];
    
    spec.k_range = K0;
    spec.phi_range = Phi;
    spec.kfactor = kfactor;
    
    disp(K0)
    
    %Solve4AllKInc(spec,AR,GridParam,ScattererType, path,filename);
    Solve4AllPhiInc(spec,ScattererType,path, filename,rchoice);
    %     PlotM(filename);
    
    %     x=linspace(-1,1,1000);
    %     y=x;
    %      GenerateSARImage(x,y,filename)
    %   load(filename);
    
end

function Solve4AllKInc(spec,AR,GridParam,ScattererType,path,filename)
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
    DiffOp = @Tools.DifferentialOps.HelmholtzOp;
    
    shift=[0,0];
    
    if  strcmpi(ScattererType,'shifted')
        ScattererType = 'StarShapedScatterer';
        
        a=0.9;
        b=a;
        
        r0=0.4;
        r1=1.4;
        
        %interior problem in square
        x1=-1.4;xn=1.4;
        y1=-1.4;yn=1.4;
        
        %FocalDistance = sqrt(a^2-b^2);
        %Eta0 = acosh(a/FocalDistance);
        
        shift=[1/3,1/3];
        
        Parameterization  = Tools.Parameterizations.ParametricEllipse2(struct('a',a,'b',b,'xcenter',1/3,'ycenter',1/3,'rotation',0));
        
        
        ScattererHandle  = @Tools.Scatterer.StarShapedScatterer;
        UParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization);
        
        ScattererParams  = UParams;
        ScattererParams.Stencil=9;
        
        Extension =  @Tools.Extensions.TwoTupleExtension;
        
    elseif  strcmpi(ScattererType,'ellipse')
        ScattererType = 'StarShapedScatterer';
        
        a=1;%b=a/2;
        b=a/AR;
        %exterior problem in ring
        %r0 = 0.7*b;
        %r1 = 1.8*a; %ABC set on that circle
        
        r0=0.8*b;
        r1=1.2*a;
        
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
    
    str = sprintf('\n ParSAR:Solve4AllKInc, Ellipse a=%4.2f, b=%4.2f, dx=%4.2f, dy=%4.2f  , filename=%s\n',a,b,shift(1),shift(2),[path filesep filename]);
    disp(str)
    parfor indx=1:numel(spec.k_range)
        k0=spec.k_range(indx);
        k1=spec.kfactor*k0;
        
        
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
            'Extension', Extension, 'ExtensionParams',[] , ...
            'DiffOp', DiffOp, 'DiffOpParams', [] ));
        
        s1=toc(t1);
        tmp=[ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ];%Q1,2 are only evaluated by request...
        s2=toc(t1);
        fprintf('tmp took %g\n',s2-s1)
        
        fprintf('done computing k0=%5.3f,k1=%5.3f after %6.4f secs\n',k0,k1,toc(t1));
        parsave([path filesep filename 'k0_' num2str(k0) '.k1_' num2str(k1) '.mat'], PlrGrid, CrtsGrid ,ExtPrb, IntPrb,UParams, k0,k1,a,b,shift,Basis);
    end
    
    fprintf('finished all after %6.4f secs in total\n',toc(t1));
       
end

function Solve4AllPhiInc(spec,ScattererType,path,filename,rchoice)
    t1=tic;

    ScattDataPath = [path filesep 'ScattData'];
    if 7~=exist(ScattDataPath,'dir'), mkdir(ScattDataPath); end


    parfor indx=1:numel(spec.k_range)
        k0=spec.k_range(indx);
        k1=spec.kfactor*k0;
        Nphi = numel(spec.phi_range);
        
        Kstr = ['k0_' num2str(k0) '.k1_' num2str(k1)];
          
          fileToLoad = [path filesep filename Kstr];
        [PlrGrid, CrtsGrid ,ExtPrb, IntPrb,UParams,k0,k1,a,b,shift,Basis] = parload(fileToLoad);
        
        fprintf('ParSAR:Solve4AllPhiInc, loaded Ellipse a=%4.2f, b=%4.2f, dx=%4.2f, dy=%4.2f  , filename=%s\n',a,b,shift(1),shift(2),fileToLoad);
                
        tmp=cell(1,numel(rchoice));
        [tmp{:}]=deal(zeros(1,PlrGrid.Nth));
        
        %prealocation of scattData for 'parfor', should be done in such 2 steps as of matlab R2017a (a workaround of a bug???)
         scattData =  struct('r',zeros(1,numel(rchoice)), ...
             'u', {tmp},...
             'un', {tmp}, ...
             'cn0',zeros(1,Basis.NBss0), ...
             'cn1',zeros(1,Basis.NBss1), ...
              'k0',k0, ...
              'k1',k1, ...
              'phi',0) ; 
             scattData = repmat( {scattData}, 1, Nphi);
        %scattData=cell(1,Nphi);
    %keyboard
    s1=toc(t1);
        %A=[ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ];
        [Q,R,E]=qr([ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ]);%,'matrix');
                s2=toc(t1);
                %fprintf('qr decomp. took %d\n',s2-s1)
        
        for jndx = 1:Nphi
            l1=toc(t1);
            IncAng = spec.phi_range(jndx);%*pi/180;
            
            UincParams = UParams;
            UincParams.r = ExtPrb.Scatterer.r;
            %UincParams  = UParams;struct('ScattererType','circle','r',ExtPrb.Scatterer.r);
            
            rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
            uinc = Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k0);
            rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
            
            s1=toc(t1);
            %cn = A \ rhs;
            cn = (E*(R\(Q'*rhs)));
            s2=toc(t1);
		%fprintf('qr solve took %d\n',s2-s1)
            
            if 1
                % Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                % Intxi(IntPrb.GridGamma) = [IntPrb.W{1}(IntPrb.GridGamma,:),IntPrb.W{2}(IntPrb.GridGamma,:)]*cn;
                % Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(PlrGrid.Nr,PlrGrid.Nth,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}(ExtPrb.GridGamma,:),ExtPrb.W{2}(ExtPrb.GridGamma,:)]*cn ;%- uinc;
                
                UincParams = UParams;
                UincParams.r = PlrGrid.R;
                %UincParams  = struct('ScattererType','circle','r',PlrGrid.R);
                uinc = Uinc(UincParams,PlrGrid.Theta,IncAng,k0);
                
                Extu = ExtPrb.P_Omega(Extxi,uinc ) - uinc;%we want the scatterer field only here
                
                
                %create SAR image
                if  strcmpi(ScattererType,'circle')
                    s= ExtPrb.Scatterer.th;
                else
                    s= ExtPrb.Scatterer.nrml_t;
                end
                
                for mm = 1:numel(rchoice)
                    nn = numel(PlrGrid.r) - rchoice(mm);
                    
                    scattData{jndx}.r(mm) = PlrGrid.r(nn);
                    
                    scattData{jndx}.u{mm} = Extu(nn,:);
                    
                    % SECOND ORDER
                    %un = (U(n+1,:) - U(n-1,:)) / (2 * Grid.dx);
                    
                    % FOURTH ORDER
                    scattData{jndx}.un{mm}= (8 * (Extu(nn+1,:) - Extu(nn-1,:)) - (Extu(nn+2,:) - Extu(nn-2,:)) ) / (12 * PlrGrid.dx);
                end
                
            end
            
            %save it all
            scattData{jndx}.cn1 = cn( (Basis.NBss0+1):end);
            scattData{jndx}.cn0 = cn(1:Basis.NBss0);
            %scattData{jndx}.k=k0;
            %scattData{jndx}.k1=k1;
            scattData{jndx}.phi = IncAng;
        end
        
        parsave2([ScattDataPath filesep Kstr],scattData);

			
        %t3=toc(t1);
        fprintf('done computing diff inc angles for k0=%5.3f,k1=%5.3f after %6.4f secs\n',k0,k1,toc(t1));
        
        
    end
    
    fileToLoad = [path filesep filename 'k0_' num2str(spec.k_range(1)) '.k1_' num2str(spec.kfactor*spec.k_range(1))  '.mat'];
    load(fileToLoad);
    
    fileToSave = [ScattDataPath filesep filename 'k' num2str(spec.k_range(1))  '_k' num2str(spec.k_range(end)) '.mat'];
    str = sprintf('Ellipse a=%4.2f, b=%4.2f, dx=%4.2f, dy=%4.2f  , filename=%s',a,b,shift(1),shift(2),  fileToSave  );
    
    M.meta = str;
    M.Basis = Basis;
       %%M.PlrGrid=PlrGrid;
    M.spec=spec;
    M.spec.theta=PlrGrid.theta;
    
    save(fileToSave,'M','a','b','-v7.3');%, 'PlrGrid', 'CrtsGrid' ,'ExtPrb', 'IntPrb');
    
   % fprintf(' saving took %6.4f secs\n',toc(t1)-t3);
    %} 
    
    fprintf('finished all after %6.4f secs in total\n',toc(t1));
   
end

function CreateSARdata(spec,ScattererType,path,filename,rchoice)
    t1=tic;

    for indx=1:numel(spec.k_range)
        k0=spec.k_range(indx);
        k1=spec.kfactor*k0;
        
        ScattDataPath = [path filesep ScattData];
        if 7~=exist(ScattDataPath,'dir'), mkdir(ScattDataPath); end
        
        load([path filesep filename 'k0_' num2str(k0) '.k1_' num2str(k1) '.mat']);%, PlrGrid, CrtsGrid ,ExtPrb, IntPrb,UParams, k0,k1);
        str = sprintf('\n ParSAR:Solve4AllPhiInc, loaded Ellipse a=%4.2f, b=%4.2f, dx=%4.2f, dy=%4.2f  , filename=%s\n',a,b,shift(1),shift(2),filename);
        disp(str)
        
        for jndx = 1:numel(spec.phi_range)
            
            IncAng = spec.phi_range(jndx);%*pi/180;
            
            UincParams = UParams;
            UincParams.r = ExtPrb.Scatterer.r;
            %UincParams  = UParams;struct('ScattererType','circle','r',ExtPrb.Scatterer.r);
            
            rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
            uinc = Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k0);
            rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
            
            cn = [ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ] \ rhs;
            
            if 1
                % Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                % Intxi(IntPrb.GridGamma) = [IntPrb.W{1}(IntPrb.GridGamma,:),IntPrb.W{2}(IntPrb.GridGamma,:)]*cn;
                % Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(PlrGrid.Nr,PlrGrid.Nth,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}(ExtPrb.GridGamma,:),ExtPrb.W{2}(ExtPrb.GridGamma,:)]*cn ;%- uinc;
                
                UincParams = UParams;
                UincParams.r = PlrGrid.R;
                %UincParams  = struct('ScattererType','circle','r',PlrGrid.R);
                uinc = Uinc(UincParams,PlrGrid.Theta,IncAng,k0);
                
                Extu = ExtPrb.P_Omega(Extxi,uinc ) - uinc;%we want the scatterer field only here
                
                %create SAR image
                if  strcmpi(ScattererType,'circle')
                    s= ExtPrb.Scatterer.th;
                else
                    s= ExtPrb.Scatterer.nrml_t;
                end
                
                for mm = 1:numel(rchoice)
                    nn = numel(PlrGrid.r) - rchoice(mm);
                    
                    M.scattData{jndx,indx}.r(mm) = PlrGrid.r(nn);
                    
                    M.scattData{jndx,indx}.u{mm} = Extu(nn,:);
                    
                    % SECOND ORDER
                    %un = (U(n+1,:) - U(n-1,:)) / (2 * Grid.dx);
                    
                    % FOURTH ORDER
                    M.scattData{jndx,indx}.un{mm}= (8 * (Extu(nn+1,:) - Extu(nn-1,:)) - (Extu(nn+2,:) - Extu(nn-2,:)) ) / (12 * PlrGrid.dx);
                end
                
            end
            
            %save it all
            M.scattData{jndx,indx}.cn1 = cn( (Basis.NBss0+1):end);
            M.scattData{jndx,indx}.cn0 = cn(1:Basis.NBss0);
            M.scattData{jndx,indx}.k=k0;
            M.scattData{jndx,indx}.phi = IncAng;
        end
        
        t3=toc(t1);
        fprintf('done computing diff inc angles for k0=%5.3f,k1=%5.3f after %6.4f secs\n',k0,k1,toc(t1));
        
        
    end
    
    M.meta = str;
    M.Basis = Basis;
    %M.PlrGrid=PlrGrid;
    M.spec=spec;
    M.spec.theta=PlrGrid.theta;
    
    save(filename,'M','a','b','-v7.3');%, 'PlrGrid', 'CrtsGrid' ,'ExtPrb', 'IntPrb');
    
    fprintf(' saving took %6.4f secs\n',toc(t1)-t3);
    
    
    fprintf('finished all after %6.4f secs in total\n',toc(t1));
    
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

function PlotM(filename)
    %not supposed to work after all changes?....
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

function ts=timestr()
    [hr,mn,sc]=hms(datetime);
    ts = sprintf('%s%d%d%d',date,hr,mn,sc);
end

function s = firstUpper(str)
    s = regexprep(str,'(\<[a-z])','${upper($1)}');
end

function parsave2(filename,scattData)
    save([filename '.mat'],'-v7.3', 'scattData');
end

function parsave(filename,PlrGrid, CrtsGrid ,ExtPrb, IntPrb,UParams,k0,k1,a,b,shift,Basis)
    save([filename '.mat'],'-v7.3', 'PlrGrid', 'CrtsGrid' ,'ExtPrb', 'IntPrb','UParams','k0','k1','a','b','shift','Basis');
end

function [PlrGrid, CrtsGrid ,ExtPrb, IntPrb,UParams,k0,k1,a,b,shift,Basis] = parload(filename)
    load([filename '.mat'],'PlrGrid', 'CrtsGrid' ,'ExtPrb', 'IntPrb','UParams','k0','k1','a','b','shift','Basis');
end