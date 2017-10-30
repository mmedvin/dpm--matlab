function RunSoliton

R0 =1;
H0 = linspace(2,6,120);

n=4;
Solve4AllH(R0,H0,n);
FindMax(R0,H0,n);

end

function Solve4AllH(R0,H0,n)

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




k =sqrt(3);

Problem = 'Neumann';
fprintf('Grid: r0=%f, r1=%f  \n',r0,r1);


fprintf('%s problem, cmpr using grid convergence, scatterer is circle, Basis is Fourier \n',Problem);

ErrPre = 0;

for indx = 1:numel(H0)
    ExParams{indx}  = struct('ScattererType','circle','r',R0, 'H', H0(indx));
    
    f1      = @(phi) -Uinc(ExParams{indx},phi,k);
    dfdn    = @(phi) -dnUinc(ExParams{indx},phi,k);
    
    Basis{indx} = Tools.Basis.FourierBasis.BasisHelper(f1,dfdn,[45,45]);%,[1e-06,1e-06]);%105);
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
        xi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
        
        if strcmpi(Problem , 'Dirichlet')
            cn1 =( Q1 \ ( -Q0*Basis{indx}.cn0 )) ;
            xi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*Basis{indx}.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
        elseif strcmpi(Problem , 'Neumann')
            cn0 =( Q0 \ ( -Q1*Basis{indx}.cn1 )) ;
            xi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*Basis{indx}.cn1;
        else
            error('Solving only Dirichlet or Neumann problem')
        end
        
        u = ExtPrb.P_Omega(xi);
        
        
        %%%%%%%%%%%%%%
        
        
        % % % % % % % % % % % % % % % %
        % Comparison
        % % % % % % % % % % % % % % % %
        
        ExParams2 = ExParams{indx};
        ExParams2.r=Grid.R;
        
        TotalField{indx} = u + Uinc(ExParams2,Grid.Theta,k);
        
    end
    
    fprintf('\n');
    


    Mp = ExtPrb.Scatterer.Mp;
    Nm = ExtPrb.Scatterer.Nm;
    save(['Soliton' num2str(n) '.mat'],'TotalField','Grid', 'Mp','Nm');

end

function FindMax(R0,H0,n)
    
    load(['Soliton' num2str(n) '.mat'],'TotalField','Grid', 'Mp');%,'Nm');
    
    MovieName = 'SolitonMovie';
    VideoWriterObj = VideoWriter([MovieName '.avi']);
    VideoWriterObj.FrameRate = 1;
    open(VideoWriterObj);

    
    for indx = 1:numel(H0)
        
        TF = TotalField{indx};
        
        MaxVal = max( full(abs(TF(:))) );
        Loc = find(MaxVal - abs(TF(:))<10e-13 );
        MaxR = Grid.R(Loc);
        MaxTh = Grid.Theta(Loc);
        MaxX = MaxR.*cos(MaxTh);
        MaxY = MaxR.*sin(MaxTh);
        
        S.H(indx)=H0(indx);
        S.val(indx) = MaxVal(1);
       % S.r(indx) = MaxR(1);
       % S.th(indx) = MaxTh(1);
        S.x(indx) = MaxX(1);
        S.y(indx) = MaxY(1);
        
        
        if 1%want2plot
            [XExt,YExt,tExtu] = Tools.Common.MaskU(Grid,Mp,TF);
            
            str1 = sprintf('R=%.2f, H=%.2f, MaxVal=%.2f',R0,H0(indx),MaxVal);
            str2 = sprintf('movie/R=%.2f_H=%.2f',R0,H0(indx));
            
            figure(1); Tools.Common.mContour(XExt,YExt,full( abs(tExtu)),[str1 ' - total field, abs '] ); saveas(figure(1),[str2 '_abs.jpg'],'jpg')
            
            %text(MaxX,MaxY,'\diamondsuit','VerticalAlignment','middle','HorizontalAlignment','center','Color','red');
            
            frame = getframe(1);
%             img = frame2im(frame);
            
            writeVideo(VideoWriterObj,frame);
            
%             [A,map] = rgb2ind(img,256);
%             DelayTime =2;
%             
%             
%             
%             if indx == 1
%                 imwrite(A,map,MovieName,'gif','LoopCount',Inf,'DelayTime',DelayTime);
%             else
%                 imwrite(A,map,MovieName,'gif','WriteMode','append','DelayTime',DelayTime);
%             end
            
            
            %        figure(2); Tools.Common.mContour(XExt,YExt,full(real(tExtu)),[str1 ' - total field, real' ]); saveas(figure(2),[str2 '_real.jpg'],'jpg')
            %        figure(3); Tools.Common.mContour(XExt,YExt,full(imag(tExtu)),[str1 ' - total field, imag' ]); saveas(figure(3),[str2 '_imag.jpg'],'jpg')
            
        end
        
        fprintf('n=%d\t H=%d\t MaxVal=%-5.2f\n',n,H0(indx),MaxVal);
        
        for m=1:numel(MaxR)
            fprintf('Val:|(%-10.6f,%-10.6f)|=%-10.6f\t at (r, theta)=(%-5.2f,%-5.2f)=(x,y)=(%-5.2f,%-5.2f)\n',real(TF(Loc(m))),imag(TF(Loc(m))),abs(TF(Loc(m))),MaxR(m),MaxTh(m)*180/pi,MaxX(m),MaxY(m));
        end  
        fprintf('\n')
    end
    
    close(VideoWriterObj)
    
    figure(2)
    plot(S.H,S.val)
    title('max val as function of H')
    saveas(figure(2),'MaxValAsfH.jpg','jpg')
    
    figure(3)   
    plot(S.x,S.y,'*')
    for indx = 1:3:numel(H0)
        text(S.x(indx),S.y(indx),num2str(S.H(indx),3),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8)
    end
    
    title('(x(H),y(H)), the lable is H')
    
     saveas(figure(3),[str2 '.jpg'],'jpg')
        
     save(['Soliton' num2str(n) 'S.mat'],'S');%,'Nm');

end

function AsymptoticSolution(R,H)

    k=sqrt(3);
    a0=-0.5*(R/H) * besselj(1,k*H)/(k*besselh(1,2,k*R));
    b1=-0.5 * besselj(1,k*H)/(0.5*k*(besselh(0,2,k*R) - besselh(2,2,k*R)));
    a2= 0.5*(R/H) * besselj(1,k*H)/(0.5*k*(besselh(1,2,k*R) - besselh(2,2,k*R)));

    
    load(['Soliton' num2str(n) '.mat'],'Grid', 'Mp');

    r =Grid.R;
    th=Grid.Theta;
    
    ExParams2 = struct('ScattererType','circle','r',r, 'H', H);
    
    u = a0*besselh(0,2,k*r) + b1*sin(th).*besselh(1,2,k*r) + a2*cos(2*th).*besselh(2,2,k*R) + Uinc(ExParams2,Grid.Theta,k);
    
    
    
    [XExt,YExt,tExtu] = Tools.Common.MaskU(Grid,Mp,u);
    
   % str1 = sprintf('R=%.2f, H=%.2f, MaxVal=%.2f',R0,H0(indx),MaxVal);
  %  str2 = sprintf('movie/R=%.2f_H=%.2f',R0,H0(indx));
    
   % figure(1); Tools.Common.mContour(XExt,YExt,full( abs(tExtu)),[str1 ' - asymptotic total field, abs '] ); 
    %saveas(figure(1),[str2 '_abs.jpg'],'jpg')
    
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

function u = Uinc(Params,phi,k)

x = Params.r .* cos(phi);
y = Params.r .* sin(phi);
r = x.^2 + (y-Params.H).^2;

%x = sqrt(H2 + r2 - 2*Params.H*Params.r*sin(t));

u = besselj(0,k*r);

end

function un = dnUinc(Params,phi,k)

x = Params.r .* cos(phi);
y = Params.r .* sin(phi);
r = x.^2 + (y-Params.H).^2;

%x = sqrt(H2 + r2 - 2*Params.H*Params.r*sin(t));

un = k*(r - Params.H.*sin(phi)).* besselj(1,k*r)./r;

end

