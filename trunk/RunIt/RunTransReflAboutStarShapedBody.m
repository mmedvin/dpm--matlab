function RunTransReflAboutStarShapedBody
    
    want2plot = 0;

    ChebyshevRange = struct('a',-pi,'b',pi);%don't change it
    
    a=1;%2.5;
    b=a/2;%0.8;
    
    %r0 = 0.7*b;
    %r1 = 1.8*a;
    
r0=0.3;
r1=1.2;

    x1=-1.2;xn=1.2;
    %y1=-1.2;  yn=1.2;
    y1=-0.7;  yn=0.7;
   
  %  R0 =0.7;
   
    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
  
    IncAng = 0;
    IncAng = IncAng*pi/180;    
    
    %doesn't expected to work Parameterization  = Tools.Parameterizations.ParametricHeart(struct('a',13/16,'b',-5/16,'c',-2/16,'d',-1/16,'e',1,'p',3));
    Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b));
    %Parameterization  = Tools.Parameterizations.ParametricKite(struct('a',1,'b',.65*2,'c',1.5));
    %Parameterization  = Tools.Parameterizations.ParametricSubmarine(struct('a',1,'b',1/2,'c',0,'p',200));
    %Parameterization  = Tools.Parameterizations.ParametricStar();
    
    
%     Problem = 'Dirichlet'; % 'Dirichlet' or 'Neumann'   
%     fprintf('Solving %s defraction problem, comparing using grid convergance, data is PlaneWave, scatterer is circle\n',Problem);
     

dbk=dbstack();

kin = [.1,3,15,30];
kex = [.1,5, 10];
        
    for ki = 1 %1:3
        
        ErrIntPre = 0;         ErrExtPre = 0;   ErrTotPre = 0;
        
    ExtWaveNumberAddParams = struct('k',kex(ki),'r0',1.6);
    IntWaveNumberAddParams = struct('k',kin(ki),'r0',1.6);
    %k=ExtWaveNumberAddParams.k;
    
        
        %UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',Eta0, 'Vark',true); % ,false);%
        UincParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization);
            
        tst_k = kin(ki);%max(kin(ki),kex(ki));
        
        f1      = @(phi) Uinc(UincParams,phi,IncAng,tst_k);
        dfdn    = @(phi) detaUinc(UincParams,phi,IncAng,tst_k);
            

        Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn,60);
        %Basis = Tools.Basis.ChebyshevBasis.BasisHelper(f1,dfdn,ChebyshevRange);

        fprintf('%s, Grid: x1=%f, xn=%f, y1=%f, yn=%f \n %s \n kin=%d kex=%d M=%d \n', dbk(1).name,x1,xn,y1,yn, Parameterization.Print, kin(ki),kex(ki), Basis.M);
    
        nmax=3;
        for n=0:nmax %run different grids
            tic
            %build grid
            
            p=4;%6;%3;%1;
            Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
            Nx=2^(n+p)+1;	Ny=2^(n+p)+1;
           
%             BasisIndices        = -M:M;
            PlrGrid                = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
            WaveNumberHandle = @Tools.Coeffs.ConstantWaveNumber;
                      
            ScattererHandle  = @Tools.Scatterer.StarShapedScatterer;
            ScattererParams  = UincParams;
            
			CollectRhs=0;
			
            ExtPrb = Solvers.ExteriorSolver ...
                (Basis,PlrGrid,WaveNumberHandle,ExtWaveNumberAddParams,ScattererHandle,ScattererParams,CollectRhs);

            CrtsGrid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);

            %WaveNumberHandle = @Tools.Coeffs.WaveNumberElliptical;
            
            CollectRhs=1;
            IntPrb = Solvers.InteriorHomoSolver ...
                (Basis,CrtsGrid,WaveNumberHandle,IntWaveNumberAddParams,ScattererHandle,ScattererParams,CollectRhs);
         
            if 1
                IntQ = IntPrb.Q;
                ExtQ = ExtPrb.Q;%[ExtPrb.Q0,-ExtPrb.Q1]; %
                
                
                %UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',ExtPrb.Scatterer.eta, 'Vark',false);
                UincParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization,'r',ExtPrb.Scatterer.r);

                rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
                % rhs(numel(IntPrb.GridGamma)+1:end,1)= -Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
                uinc = Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,kex(ki));
                rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
                
                cn = [ IntQ ; ExtQ ] \ rhs;
                
                Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*cn;
                Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = ExtPrb.W(ExtPrb.GridGamma,:)*cn ;
                
                UincParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization,'r',PlrGrid.R);
                Uinc = Uinc(UincParams,PlrGrid.Theta,IncAng,kex(ki));
                
                Extu = ExtPrb.P_Omega(Extxi,Uinc ); 
                
            else
                
              
                cn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*Basis.cn0 )) ; 
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
                Extu = ExtPrb.P_Omega(Extxi,zeros(Nr,Nth-1));
                
                
                cn1 =( IntPrb.Q1 \ ( -IntPrb.Q0*Basis.cn0 )) ;
                
                Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Basis.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1;
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
                ErrExt =norm(tmp(:),inf);
                
                tmp = Intu(1:2:end,1:2:end)-Intu1(1:2:end,1:2:end);
                ErrInt =norm(tmp(:),inf);
                
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
                plot(a*cos(th),b*sin(th),'k','LineWidth',2)
                hold off;
                
                figure(2)
                pcolor([XInt,XExt],[YInt,YExt],real(full([tIntu,tExtu])))
                
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
                plot(a*cos(th),b*sin(th),'k','LineWidth',2)
                hold off;
                
                figure(3)
                pcolor([XInt,XExt],[YInt,YExt],imag(full([tIntu,tExtu])))
                
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
                plot(a*cos(th),b*sin(th),'k','LineWidth',2)
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


function uinc = Uinc(Params,phi,IncAng,k)  
    
    if strcmpi(Params.ScattererType,'ellipse')
         x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
         y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        x = Params.r .* cos(phi);
        y = Params.r .* sin(phi);
    elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
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
    
    if strcmpi(Params.ScattererType,'ellipse')
        dx = Params.FocalDistance * sinh(Params.eta) .* cos(phi);
        dy = Params.FocalDistance * cosh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        dx = cos(phi);
        dy = sin(phi);
    elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
        [x,dx] = Params.Parameterization.XHandle.Derivatives(phi);
        [y,dy] = Params.Parameterization.YHandle.Derivatives(phi);
    end
    
   
    
    uinc = Uinc(Params,phi,IncAng,k)  ;
    
    duinc = 1i .* k .*  uinc .* (dx.*cos(IncAng) + dy.*sin(IncAng));
 %   h = FocalDist*sqrt(sinh(eta).^2 + sin(phi).^2);
   % duinc = duinc./h;

   if strcmpi(Params.ScattererType,'StarShapedScatterer')
	   h = sqrt(dx.^2 + dy.^2);
duinc = 1i .* k .*  uinc .* (dy.*cos(IncAng) - dx.*sin(IncAng))./h;
   end
   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



