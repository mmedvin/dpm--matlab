function RunTransReflAboutStarShapedBody
    
    want2plot = 0;

    ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
    
    a=1;%2.5;
    b=a/2;%/2;%0.8;
    
    %r0 = 0.7*b;
    %r1 = 1.8*a;
    
%r0=0.8;
%r1=2.2;

	%x1=-pi; xn=pi;
	%y1=-pi; yn=pi;

    %x1=-1.2;xn=1.2;
    %y1=-1.2;  yn=1.2;
    %y1=-0.7;  yn=0.7;
   

    % Kite
%    x1=-1.7;xn=1.2;
%    y1=-1.7;  yn=1.7;
%    r0=0.8;
%    r1=2.2;

   % star
    x1=-1.7; xn=1.7;
    y1=-1.7; yn=1.7;
    r0=0.3;
    r1=2.2;


   % submarine
%     x1=-2.2; xn=2.2;
%     y1=-0.6; yn=1.2;
%     r0=0.3;
%     r1=2.2;

   
  %  R0 =0.7;
   
    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
  
    IncAngD = 0;
    IncAng = IncAngD*pi/180;    
    
    %doesn't expected to work Parameterization  = Tools.Parameterizations.ParametricHeart(struct('a',13/16,'b',-5/16,'c',-2/16,'d',-1/16,'e',1,'p',3));
    %Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b));
    %Parameterization  = Tools.Parameterizations.ParametricKite(struct('a',1,'b',.65*2,'c',1.5));
    %Parameterization  = Tools.Parameterizations.ParametricSubmarine(struct('a',1,'b',1/2,'c',0,'p',150));
    %Parameterization  = Tools.Parameterizations.ParametricSubmarine(struct('a',1.8,'b',1.8/5,'c',1,'p',150));
   % Parameterization  = Tools.Parameterizations.ParametricStar();
   
     Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b,'xcenter',.05,'ycenter',0.05,'rotation',pi/2));

    
%     Problem = 'Dirichlet'; % 'Dirichlet' or 'Neumann'   
%     fprintf('Solving %s defraction problem, comparing using grid convergance, data is PlaneWave, scatterer is circle\n',Problem);
     

dbk=dbstack();

kin = [3 ,10, 15,  20];
kex = [1 ,5 ,  5,  10];
        
    for ki = 1 %1:3
        
        ErrIntPre = 0;         ErrExtPre = 0;   ErrTotPre = 0;
                    
        %UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',Eta0, 'Vark',true); % ,false);%
        UincParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization);
            
        tst_k = max(kin(ki),kex(ki));
  	%tst_k = kin(ki) + kex(ki);

        f1      = @(phi) Uinc(UincParams,phi,IncAng,tst_k);
        dfdn    = @(phi) detaUinc(UincParams,phi,IncAng,tst_k);
            

        Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn,[1e-5,1e-5]);
	%Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn,fix(Basis.M/2));
        %Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f1,dfdn,ChebyshevRange);

        fprintf('%s, IncAngD: %f, Grid:  x1=%f, xn=%f, y1=%f, yn=%f, r0=%f, r1=%f \n %s \n kin=%d kex=%d NBss0=%d NBss1=%d \n', dbk(1).name,IncAngD,x1,xn,y1,yn,r0,r1, Parameterization.Print, kin(ki),kex(ki), Basis.NBss0, Basis.NBss1);
    
        nmax=4;%3;
        for n=1:nmax %run different grids
            tic
            %build grid
            
            p=5;%4;%6;%3;%1;
            Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
            Nx=2^(n+p)+1;	Ny=2^(n+p)+1;
           
%             BasisIndices        = -M:M;
            PlrGrid                = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
                      
            ScattererHandle  = @Tools.Scatterer.StarShapedScatterer;
            ScattererParams  = UincParams;
            ScattererParams.Stencil=9;
            
            ExtPrb = Solvers.ExteriorSolver( struct(...
                     'Basis',Basis, ...
                     'Grid', PlrGrid, ...
                     'CoeffsHandle', @Tools.Coeffs.ConstantWaveNumber, ... 
                     'CoeffsParams',  struct('k',kex(ki),'r0',1.6), ...
                     'ScattererHandle',ScattererHandle, ...
                     'ScattererParams', ScattererParams, ...
                     'CollectRhs',0, ... %i.e. no 
                     'Extension', @Tools.Extensions.TwoTupleExtension, ...
                     'ExtensionParams',[] ...
                     ));

            CrtsGrid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
            
            IntPrb = Solvers.InteriorHomoSolver( struct(...
                     'Basis',Basis, ...
                     'Grid', CrtsGrid, ...
                     'CoeffsHandle', @Tools.Coeffs.WaveNumberStarShaped, ... %ConstantWaveNumber;
                     'CoeffsParams', struct('k',kin(ki),'r0',1.6), ...
                     'ScattererHandle',ScattererHandle, ...
                     'ScattererParams', ScattererParams, ...
                     'CollectRhs',1, ... %i.e. yes
                     'Extension', @Tools.Extensions.TwoTupleExtension, ...
                     'ExtensionParams',[] , ...
                     'DiffOp'            , @Tools.DifferentialOps.HelmholtzOp, ...
                     'DiffOpParams'      , [] ...
                     )); ...
         
            if 1
                
                
                %UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',ExtPrb.Scatterer.eta, 'Vark',false);
                UincParams1 = struct('ScattererType','StarShapedScatterer','r',ExtPrb.Scatterer.r);%,'Parameterization',Parameterization

                rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
                % rhs(numel(IntPrb.GridGamma)+1:end,1)= -Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
                
                uinc = Uinc(UincParams1,ExtPrb.Scatterer.th,IncAng,kex(ki));
                rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
                
                %[uinc,uinc_t,uinc_tt,uinc_3t,uinc_4t] = Uinc(    UincParams,ExtPrb.Scatterer.BasisArg,IncAng,kex(ki));
                %[duinc,duinc_t,duinc_tt]              = detaUinc(UincParams,ExtPrb.Scatterer.BasisArg,IncAng,kex(ki));
                 
                    %Xi0 = Tools.Basis.BasisFunctionWD(); 
                    %Xi1 = Tools.Basis.BasisFunctionWD();
                    
                    %Xi0.xi0        = uinc;
                    %Xi0.xi0t       = uinc_t;
                    %Xi0.xi0tt      = uinc_tt;
                    %Xi0.xi0ttt     = uinc_3t;
                    %Xi0.xi0tttt    = uinc_4t;
                    
                    %Xi1.xi0        = duinc;
                    %Xi1.xi0t       = duinc_t;
                    %Xi1.xi0tt      = duinc_tt;
                    %Xi1.xi0ttt     = 0;
                    %Xi1.xi0tttt    = 0;
                
                %rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(Xi0,Xi1);

                cn = [ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ] \ rhs;
                
                Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = [IntPrb.W{1}(IntPrb.GridGamma,:),IntPrb.W{2}(IntPrb.GridGamma,:)]*cn;
                Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}(ExtPrb.GridGamma,:),ExtPrb.W{2}(ExtPrb.GridGamma,:)]*cn ;%- uinc;
                
                 UincParams2 = struct('ScattererType','StarShapedScatterer','r',PlrGrid.R);%,'Parameterization',Parameterization
                Uinc = Uinc(UincParams2,PlrGrid.Theta,IncAng,kex(ki));
                
                Extu = ExtPrb.P_Omega(Extxi,Uinc ); 
                
            else
                
                cn0 = Basis.cn0;
                cn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*Basis.cn0 )) ; 
		%cn0 =( ExtPrb.Q0 \ ( -ExtPrb.Q1*Basis.cn1 )) ; 
		%cn1 = Basis.cn1;
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
                Extu = ExtPrb.P_Omega(Extxi,zeros(Nr,Nth-1));
                
                cn0 = Basis.cn0;
                cn1 =( IntPrb.Q1 \ ( -IntPrb.Q0*Basis.cn0 )) ;
		%cn0 =( IntPrb.Q0 \ ( -IntPrb.Q1*Basis.cn1 )) ; 
                %cn1 = Basis.cn1;
                
                Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1;
                Intu = IntPrb.P_Omega(Intxi);
                
                
%                 Intu = zeros(Nx,Ny); %IntPrb.P_Omega(Intxi);
                
                
            end
            
            t=toc;
            
            %%%%%%%%%%%%%%
            
            % % % % % % % % % % % % % % % %
            % Comparison
            % % % % % % % % % % % % % % % %
            if n > 1
                Extu1= spalloc(Nr,Nth-1,nnz(Extu0));
                Extu1(ExtPrb.Scatterer.Nm) = Extu0(ExtPrb.Scatterer.Nm);
                
                Intu1= spalloc(Nx,Ny,nnz(Intu0));
                Intu1(IntPrb.Np) = Intu0(IntPrb.Np);
                
                
                tmp = Extu(1:2:end,1:2:end)-Extu1(1:2:end,1:2:end);
tmp2 = 1;%Extu(1:2:end,1:2:end);
                ErrExt =norm(tmp(:),inf)/norm(tmp2(:),inf);;
                
                tmp = Intu(1:2:end,1:2:end)-Intu1(1:2:end,1:2:end);
tmp2 = 1;%Intu(1:2:end,1:2:end);
                ErrInt =norm(tmp(:),inf)/norm(tmp2(:),inf);
                
                ErrTot = max(ErrInt,ErrExt);
                
               % fprintf('kex=%d,kin=%d,M=%d,Nplr=%-5dx%d\t, Ncrt=%-5dx%d\t ExtErr=%d\t IntErr=%d\t time=%d\n',k,k+dk,Basis.M, Nr,Nth,Nx,Ny,full(ExtErr(n)),full(IntErr(n)),t);
               %kex=%d,kin=%d,M=%d,
               %kex(ki),kin(ki),Basis.M,
               fprintf('N=%-10dx%-10d\t ErrExt=%d rate=%-6.2f  ErrInt=%d rate=%-6.2f  ErrTot=%d rate=%-6.2f\t time=%d\n',...
                         Nr,Nth,ErrExt,log2(ErrExtPre/ErrExt),ErrInt,log2(ErrIntPre/ErrInt),ErrTot,log2(ErrTotPre/ErrTot),t);

       		ErrIntPre = ErrInt ; 
        	ErrExtPre = ErrExt ;
            ErrTotPre = ErrTot;

            end
            
            Extu0=spalloc(Nr*2-1,Nth*2-2,nnz(Extu));
            Extu0(1:2:end,1:2:end)=Extu;
            
            Intu0=spalloc(Nx*2-1,Ny*2-1,nnz(Intu));
            Intu0(1:2:end,1:2:end)=Intu;

            if want2plot && n>0 %|| n==nmax

               % XExt = PlrGrid.R .* cos(PlrGrid.Theta);
               % YExt = PlrGrid.R .* sin(PlrGrid.Theta);
               % XExt(ExtPrb.Scatterer.Mp)=NaN;
               % YExt(ExtPrb.Scatterer.Mp)=NaN;
                
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
               tExtu(Nm)=Extu(Nm);
               tExtu(:,Nth)=tExtu(:,1);
                
               Np=IntPrb.Scatterer.Np;
               tIntu = ones(size(CrtsGrid.X)).*NaN;
                tIntu(Np) = Intu(Np);
                %tIntu(IntPrb.Scatterer.Mm)=0;%NaN;
                
                XInt = ones(size(CrtsGrid.X)).*NaN;
                YInt = ones(size(CrtsGrid.X)).*NaN;
                
                XInt(Np) = CrtsGrid.X(Np);
                YInt(Np) = CrtsGrid.Y(Np);
                % XInt(IntPrb.Scatterer.Mm)=NaN;
                %YInt(IntPrb.Scatterer.Mm)=NaN;
                
                %filename = sprintf('TFElpsAR%dRing%d%dkin%dkex%dgrd%d',fix(10*a/b),fix(10*r0),fix(10*r1),IntWaveNumberAddParams.k,ExtWaveNumberAddParams.k,Nr);
                %ExScat = ExtPrb.Scatterer;
                %IntScat = IntPrb.Scatterer;
                %save(filename,'filename','a','b','XInt','XExt','YInt','YExt','tIntu','tExtu','Extu','Intu','Nm','Np','PlrGrid','CrtsGrid','ExScat','IntScat');
                
               % continue
                
                th=0:0.0001:2*pi;
                
                figure(1)
                pcolor([XInt,XExt],[YInt,YExt],abs(full([tIntu,tExtu])))
                
                colormap jet
                title('total field, abs');
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
                pcolor([XInt,XExt],[YInt,YExt],real(full([tIntu,tExtu])))
                
                colormap jet
                title('total field, real');
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
                pcolor([XInt,XExt],[YInt,YExt],imag(full([tIntu,tExtu])))
                
                colormap jet
                title('total field, imag');
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
                
              fprintf('press key');
              soundsc(sin(2000*(0:0.001:2*pi)).^24);
              pause
              fprintf('.');
            end
            
            
        end
        
        fprintf('\n');
        
       % ExtL=log2(ExtErr(1:end-1)./ExtErr(2:end))
       % IntL=log2(IntErr(1:end-1)./IntErr(2:end))
    end
    
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

