function RunSoliton(type)

R0 =1;
n=4;

%AsymptoticSolution(R0,6,n);

if type == 1

    path = [pwd '\Results\1NoTime'];
    diary([path filesep 'result.txt']);
    
    H0  = linspace(2,6,120);
    Tau0= 0;
        
    k =sqrt(3);
    
    %Solve4AllH(k,R0,H0,Tau0,n,@Uinc1,@dnUinc1, path);    

    FindMax(k,R0,H0, Tau0,n,@Uinc1,path);

    diary off
elseif type==2
    path = [pwd '\Results\1Time'];
    diary([path filesep 'result.txt']);
    
    H0 = 2.5;
    Tau0= linspace(0,5,100);
    AddParams = struct('L',2.5,'Beta',0.1);
    AddParams.Gamma = (1-AddParams.Beta^2)^(-1/2);   
    
    k = sqrt(4*AddParams.Gamma^2 - 1);
    
    Solve4AllH(k,R0,H0,Tau0,n,@Uinc2,@dnUinc2,path,AddParams);

    FindMax(k,R0,H0, Tau0,n,@Uinc2,path);
    diary off
end


end

function Solve4AllH(k,R0,H0,Tau0,n,Uinc,dnUinc,path, AddParams)

    LastKnownVersionOfTools = '2.0.0.0';
    try
        ver = Tools.Version();
        if ~strcmp(ver,LastKnownVersionOfTools)
            error('MDP:wrong version of Tools, expected version %s, found version %s',LastKnownVersionOfTools,ver);
        end
    catch err
        if strcmp(err.identifier, 'MATLAB:undefinedVarOrClass')
            error('MDP: please add parent folder to the path');
        else%if strcmp(err.identifier,'MDP:wrong version of Tools')
            sprintf(err.message);
            rethrow(err);
        end
    end
    
    % want2plot = 0;
    
    
    r0=0.8; r1=7;
    NHR=1.6;
    
    
    Problem = 'Neumann';
    fprintf('Grid: r0=%f, r1=%f  \n',r0,r1);
    
    
    fprintf('%s problem, cmpr using grid convergence, scatterer is circle, Basis is Fourier \n',Problem);
    
    if nargin == 9
        fprintf('L=%f, Beta=%f,Gamma=%f\n',AddParams.L, AddParams.Beta,AddParams.Gamma);
    end
    
    ErrPre = 0;
    
    for indx = 1:numel(H0)
        for jndx = 1:numel(Tau0)
            ExParams{indx,jndx}  = struct('ScattererType','circle','r',R0, 'H', H0(indx));
            
            if nargin == 9
                ExParams{indx,jndx}.L = AddParams.L;
                ExParams{indx,jndx}.Beta = AddParams.Beta;
                ExParams{indx,jndx}.Gamma = AddParams.Gamma;
                ExParams{indx,jndx}.Tau = Tau0(jndx);                               
            end
            
            f1      = @(phi) -Uinc(ExParams{indx,jndx},phi,k);
            dfdn    = @(phi) -dnUinc(ExParams{indx,jndx},phi,k);
            
            Basis{indx,jndx} = Tools.Basis.FourierBasis.BasisHelper(f1,dfdn,[30,50]);%[45,45]);%,[1e-06,1e-06]);%105);
        end
    end


    %build grid
    
    p=6;
    Nr=2^(n+p)+1;
    Nth=2^(n+p)+1;
    Grid = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
    
    ExtPrb =  Solvers.ExteriorSolver( struct(...
        'Basis'           , Basis{1}, ...
        'Grid'            , Grid, ...
        'CoeffsHandle'    , @Tools.Coeffs.ConstantWaveNumber, ...
        'CoeffsParams'    , struct('k',k,'r0',NHR), ...
        'ScattererHandle' , @Tools.Scatterer.PolarScatterer, ...
        'ScattererParams' , struct('r0',R0,'ExpansionType',15,'Stencil',9), ...
        'CollectRhs'      , 0, ... %i.e. not
        'Extension'       , @Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension, ...
        'ExtensionParams' , [] ...
        ));
    
    Q0 = ExtPrb.Q0;
    Q1 = ExtPrb.Q1;
    
    for indx = 1:numel(H0)
        for jndx = 1:numel(Tau0)
            xi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
            
            if strcmpi(Problem , 'Dirichlet')
                cn1 =( Q1 \ ( -Q0*Basis{indx,jndx}.cn0 )) ;
                xi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*Basis{indx,jndx}.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
            elseif strcmpi(Problem , 'Neumann')
                cn0 =( Q0 \ ( -Q1*Basis{indx,jndx}.cn1 )) ;
                xi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*Basis{indx,jndx}.cn1;
            else
                error('Solving only Dirichlet or Neumann problem')
            end
            
            ScatteredField{indx,jndx} = full(ExtPrb.P_Omega(xi));
                                   
         end
    end
    
    fprintf('\n');
    


    Mp = ExtPrb.Scatterer.Mp;
    Nm = ExtPrb.Scatterer.Nm;
    NBss0 = Basis{indx,jndx}.NBss0;
    NBss1 = Basis{indx,jndx}.NBss1;
    save([path filesep 'Soliton' num2str(n) '.mat'],'ScatteredField','Grid','ExParams', 'Mp','Nm','NBss0','NBss1');

end

function FindMax(k,R0,H0,Tau0,n,Uinc, path)
    %if nargin ==4, path = '.'; end
    
    LoadThisFileName = [path filesep 'Soliton' num2str(n) '.mat'];
    %MovieDirName = [path filesep 'movie'];
    %if ~exist(MovieDirName,'dir'), mkdir(MovieDirName), end
    
    load(LoadThisFileName);%,'ScatteredField','Grid', 'Mp','NBss0','NBss1');%,'Nm');
    
    MovieName = 'SolitonMovie';
    VideoWriterObj = VideoWriter([path filesep MovieName '.avi']);
    VideoWriterObj.FrameRate = 1;
    open(VideoWriterObj);

    
    for indx = 1:numel(H0)
        for jndx = 1:numel(Tau0)
            
            SF = ScatteredField{indx,jndx};          
            ExParams2 = ExParams{indx,jndx};
            ExParams2.r=Grid.R;
            
            IncidentField = Uinc(ExParams2,Grid.Theta,k);
            TotalField = SF + IncidentField;           
            
            S.H(indx)        = H0(indx);
            S.Tau(indx,jndx) = Tau0(jndx);
            [S.val(indx,jndx),Loc,S.x(indx,jndx),S.y(indx,jndx),MaxR,MaxTh] = Max(TotalField,Grid);
            
            [S.ival(indx,jndx),iLoc,S.ix(indx,jndx),S.iy(indx,jndx)] = Max(IncidentField,Grid);
                        
            if 1%want2plot
                [Xtot ,Ytot ,tT] = Tools.Common.MaskU(Grid,Mp,TotalField);
                [Xinc ,Yinc ,tI] = Tools.Common.MaskU(Grid,Mp,IncidentField);
                [Xscat,Yscat,tS] = Tools.Common.MaskU(Grid,Mp,SF);
                
                PlotTitle = sprintf('R=%.2f, H=%.2f, Tau=%.2f, TotMaxVal=%.2f, IncMaxVal=%.2f',R0,H0(indx),Tau0(jndx),S.val(indx,jndx),S.ival(indx,jndx));
                PlotText = sprintf('TotMaxVal=%.2f at (%-5.2f,%-5.2f),\n IncMaxVal=%.2f at (%-5.2f,%-5.2f)',...
                    S.val(indx,jndx), S.x(indx,jndx),S.y(indx,jndx), ...
                    S.ival(indx,jndx),S.ix(indx,jndx),S.iy(indx,jndx));
                
                FileName = sprintf('R=%.2f_H=%.2f_Tau=%.2f',R0,H0(indx),Tau0(jndx));
                
                h=figure(1); 
                title(PlotTitle);
                
                set(h,'units','normalized','outerposition',[0 0 1 1])

                subplot(2,2,1);plot(1);
                title(PlotTitle);
                axis([0,1,0,1]);
                text(0.5,0.5,PlotText);
                axis off;
                
                subplot(2,2,3);
                Tools.Common.mContour(Xtot ,Ytot,full( abs(tT)), '|total field|' ); 

                subplot(2,2,2);
                Tools.Common.mContour(Xinc ,Yinc,full( abs(tI)), '|incident field|' ); 
                
                subplot(2,2,4);
                Tools.Common.mContour(Xscat,Yscat,full( abs(tS)), '|scattering field|' ); 

                %saveas(figure(1),[MovieDirName filesep FileName '_abs.jpg'],'jpg')
                              
                frame = getframe(1);
                writeVideo(VideoWriterObj,frame);
                                
            end
            
            fprintf('k=%-5.2f,n=%d\t H=%d\t Tau=%d\t MaxVal=%-5.2f, NBss0=%d,NBss1=%d\n',k,n,H0(indx),Tau0(jndx),S.val(indx,jndx),NBss0,NBss1);
            

            fprintf('Val:|(%-10.6f,%-10.6f)|=%-10.6f\t at (r, theta)=(%-5.2f,%-5.2f)=(x,y)=(%-5.2f,%-5.2f)\n',...
                            real(TotalField(Loc)),imag(TotalField(Loc)),abs(TotalField(Loc)),MaxR,MaxTh*180/pi,S.x(indx,jndx),S.y(indx,jndx));
            fprintf('\n');
            
        end
    end
    
    close(VideoWriterObj)
    
    if numel(Tau0)==1
        figure(2)
        plot(S.H,S.val)
        title('max val as function of H')
        saveas(figure(2),[path filesep 'MaxValAsfH.jpg'],'jpg')
        
        figure(3)
        plot(S.x,S.y,'*')
        for indx = 1:3:numel(H0)
                text(S.x(indx),S.y(indx),num2str(S.H(indx),3),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8)
        end
        
        title('(x(H),y(H)), the label is H')
        
        saveas(figure(3),[path filesep 'MaxLocation.jpg'],'jpg')

        figure(4)
        plot(S.H,S.y,'r',S.H,S.iy,'b'); legend('scat','inc');
        title('y_{max}(H)')
        xlabel('H');
        ylabel('y');
        saveas(figure(4),[path filesep 'ymax.jpg'],'jpg')

    
    elseif numel(H0)==1
        figure(2)
        plot(S.Tau,S.val)
        title('max val as function of Tau')
        saveas(figure(2),[path filesep 'MaxValAsfTau.jpg'],'jpg')
        
        figure(3)
        plot(S.x,S.y,'*')
        for indx = 1:3:numel(Tau0)
            text(S.x(indx),S.y(indx),num2str(S.Tau(indx),3),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8)
        end
        
        title('(x(\tau),y(\tau)), the label is \tau')
        
        saveas(figure(3),[path filesep 'MaxLocation.jpg'],'jpg')
        
    end
     save([path filesep 'Soliton' num2str(n) 'S.mat'],'S');%,'Nm');

end

function AsymptoticSolution(R,H,n)

    k=sqrt(3);
    a0=-0.5*(R/H) * besselj(1,k*H)/(k*besselh(1,2,k*R));
    b1=-0.5 * besselj(1,k*H)/(0.5*k*(besselh(0,2,k*R) - besselh(2,2,k*R)));
    a2= 0.5*(R/H) * besselj(1,k*H)/(0.5*k*(besselh(1,2,k*R) - besselh(2,2,k*R)));

    
    load(['Soliton' num2str(n) '.mat'],'Grid', 'Mp');

    r =Grid.R;
    th=Grid.Theta;
    
    ExParams2 = struct('ScattererType','circle','r',r, 'H', H);
    
    u = a0*besselh(0,2,k*r) + b1*sin(th).*besselh(1,2,k*r) + a2*cos(2*th).*besselh(2,2,k*R) + Uinc1(ExParams2,Grid.Theta,k);
    
    MaxVal = max( full(abs(u(:))) );
    
    [XExt,YExt,tExtu] = Tools.Common.MaskU(Grid,Mp,u);
    
    str1 = sprintf('R=%.2f, H=%.2f, MaxVal=%.2f',R,H,MaxVal);
    str2 = sprintf('Asymptotic_R=%.2f_H=%.2f',R,H);
    
    figure(1); 
   Tools.Common.mContour(XExt,YExt,full( abs(tExtu)),[str1 ' - asymptotic total field, abs '] ); 
    saveas(figure(1),[str2 '_abs.jpg'],'jpg')
    
end

function Convergance(indx)
load('Soliton.mat')

%         MaxVal = max( full(abs(TotalField(:))) );
%         Loc = find(MaxVal - abs(TotalField(:))<10e-13 );
%         MaxR = Grid.R(Loc);
%         MaxTh = Grid.Theta(Loc);
%         MaxX = MaxR.*cos(MaxTh);
%         MaxY = MaxR.*sin(MaxTh);
%         
%         if want2plot
%             [XExt,YExt,tExtu] = Tools.Common.MaskU(Grid,ExtPrb.Scatterer.Mp,TotalField);
%             
%             str1 = sprintf('R=%.2f, H=%.2f, MaxVal=%.2f',R0,H0,MaxVal);
%             str2 = sprintf('R=%.2f_H=%.2f',R0,H0);
%             
%             figure(1); Tools.Common.mContour(XExt,YExt,full( abs(tExtu)),[str1 ' - total field, abs '] ); saveas(figure(1),[str2 '_abs.jpg'],'jpg')
%             
%             text(MaxX,MaxY,'\diamondsuit','HorizontalAlignment','center','Color','red');
%             
%             %        figure(2); Tools.Common.mContour(XExt,YExt,full(real(tExtu)),[str1 ' - total field, real' ]); saveas(figure(2),[str2 '_real.jpg'],'jpg')
%             %        figure(3); Tools.Common.mContour(XExt,YExt,full(imag(tExtu)),[str1 ' - total field, imag' ]); saveas(figure(3),[str2 '_imag.jpg'],'jpg')
%             
%         end
%         
%         if n > 1
%             u1= spalloc(Nr,Nth-1,nnz(u0));
%             u1(ExtPrb.Nm) = u0(ExtPrb.Nm);
%             
%             tmp = u(1:2:end,1:2:end)-u1(1:2:end,1:2:end);
%             
%             %etinf(n) =norm(tmp(:),inf);
%             ErrTot = norm(tmp(:),inf);
%             fprintf('k=%d,NBss0=%d,NBss1=%d,N=%-10dx%-10d\t ErrTot=%d\t rate=%-5.2f\t MaxVal=%-5.2f\t MaxX=%-5.2f\t MaxY=%-5.2f\t time=%d\n',...
%                 k,Basis.NBss0,Basis.NBss1, Nr,Nth,ErrTot,log2(ErrPre/ErrTot),MaxVal,MaxX,MaxY,t);
%             
%             ErrPre = ErrTot;
%         else
%             fprintf('k=%d,NBss0=%d,NBss1=%d,N=%-10dx%-10d\t ErrTot=NA\t\t\t\t rate=NA\t MaxVal=%-5.2f\t MaxX=%-5.2f\t MaxY=%-5.2f\t time=%d\n',k,Basis.NBss0,Basis.NBss1, Nr,Nth,MaxVal,MaxX,MaxY,t);
%         end
%         
%         u0=spalloc(Nr*2-1,Nth*2-2,nnz(u));
%         u0(1:2:end,1:2:end)=u;

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MaxVal,MaxLoc,MaxX,MaxY,MaxR,MaxTh] = Max(U,Grid)
    
    [MaxVal,MaxLoc] = max( full(abs(U(:))) );
    %Loc = find(MaxVal - abs(U(:))<10e-13 );
    MaxR = Grid.R(MaxLoc);
    MaxTh = Grid.Theta(MaxLoc);
    MaxX = MaxR.*cos(MaxTh);
    MaxY = MaxR.*sin(MaxTh);        
end



function u = Uinc1(Params,phi,k)

x = Params.r .* cos(phi);
y = Params.r .* sin(phi);
r = sqrt(x.^2 + (y-Params.H).^2);

u = besselj(0,k*r);

end

function un = dnUinc1(Params,phi,k)

x = Params.r .* cos(phi);
y = Params.r .* sin(phi);
r = sqrt(x.^2 + (y-Params.H).^2);

un = -k*(r - Params.H.*sin(phi)).* besselj(1,k*r)./r;

end



function u = Uinc2(Params,phi,k)

x = Params.r .* cos(phi);
y = Params.r .* sin(phi);
r = sqrt( (Params.Gamma^2)*(x + Params.L - Params.Tau).^2 + (y-Params.H).^2 );

u = exp(2i*Params.Beta*Params.Gamma*(x+Params.L)).*besselj(0,k*r);

end

function un = dnUinc2(Params,phi,k)

x = Params.r .* cos(phi);
y = Params.r .* sin(phi);
r = sqrt( (Params.Gamma^2)*(x + Params.L - Params.Tau).^2 + (y-Params.H).^2 );

dr = ( (cos(phi).*Params.Gamma^2 ) .* (x + Params.L - Params.Tau) + (y-Params.H).*sin(phi) ) ./ r;

un = exp(2i*Params.Beta*Params.Gamma*(x+Params.L)).* (...
    2i*Params.Beta*Params.Gamma*cos(phi).* besselj(0,k*r) ...
    - k*dr.* besselj(1,k*r) ...
    );
end

