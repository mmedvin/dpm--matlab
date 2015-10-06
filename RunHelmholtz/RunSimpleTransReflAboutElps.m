function RunSimpleTransReflAboutElps
    
    want2plot = 0;

    
    
    a=1.8;%1;%2.5;
    %b=0.8;
    
    r0 = 0.1;%0.7*b;
    r1 = 2.5; %1.8*a;
    
    
%     r0 = 0.6;
%     r1 = 2;
    
%    x1=-1.2;xn=1.2;
	x1=-2;xn=2;
%	y1=-1.2;  yn=1.2;
   
  %  R0 =0.7;
NHR = 1.6;
   for b = 0.9 % [0.18, 0.35, 0.6,0.9] %[0.9,0.6,0.35,0.18,0.15]

IntErr=[];ExtErr=[];

y1=-(b+0.2);  yn=b+0.2;

    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
  
    IncAng = 40;
    IncAng = IncAng*pi/180;    
    
    dk = 30;
    
%     Problem = 'Dirichlet'; % 'Dirichlet' or 'Neumann'   
%     fprintf('Solving %s defraction problem, comparing using grid convergance, data is PlaneWave, scatterer is circle\n',Problem);
fprintf('Trans/Refl problem about ellipse of FD=%d, ,Eta0=%d, a=%d, b=%d, AR=%d ,r0=%d,r1=%d,x=+/-%d,y+/-=%d\n',FocalDistance , Eta0 , a,b,a/b,r0,r1,xn,yn);
     


for k = 1%[1, 5,15];%15%[1,3,5,10]%[1,5,10,15,20,25]
    ExtWaveNumberParams = struct('k',k,'r0',NHR);
    IntWaveNumberParams =  struct('k',3*k,'r0',NHR);
    kmax = max(ExtWaveNumberParams.k ,IntWaveNumberParams.k );
    UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',Eta0);
    f1      = @(phi) Uinc(UincParams,phi,IncAng,kmax);
    dfdn    = @(phi) detaUinc(UincParams,phi,IncAng,kmax);
    
    Basis = Tools.Basis.FourierBasis.BasisHelper(f1,dfdn);
        
      %  [cn0ex,cn1ex,M] = FourierCoeff(f1,dfdn);
    
        nmax=4;%7;%2
        for n=0:nmax %run different grids
            tic
            %build grid
            
            p=5;
            Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
            Nx=2^(n+p)+1;	Ny=2^(n+p)+1;
           
            %BasisIndices        = -M:M;
            PlrGrid                = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
            WaveNumberHandle = @Tools.Coeffs.ConstantWaveNumber;
                                   
            ScattererHandle  = @Tools.Scatterer.EllipticScatterer;
            ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'Stencil',9);
            			
            ExtPrb = Solvers.ExteriorSolver( struct(...
                     'Basis',Basis, ...
                     'Grid', PlrGrid, ...
                     'CoeffsHandle', WaveNumberHandle, ... 
                     'CoeffsParams', ExtWaveNumberParams, ...
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
                     'CoeffsHandle', WaveNumberHandle, ... 
                     'CoeffsParams', IntWaveNumberParams, ...
                     'ScattererHandle',ScattererHandle, ...
                     'ScattererParams', ScattererParams, ...
                     'CollectRhs',1, ... %i.e. yes
                     'Extension', @Tools.Extensions.TwoTupleExtension, ...
                     'ExtensionParams',[] ...
                     ));
         
            if 1             
                UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',ExtPrb.Scatterer.eta);
                
                rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
                % rhs(numel(IntPrb.GridGamma)+1:end,1)= -Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
                uinc = Uinc(UincParams,ExtPrb.Scatterer.phi,IncAng,k);
                rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
                
                cn = [ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ] \ rhs;
                
                %             cn0 = cn(1:2*M+1);
                %             cn1 = cn(2*M+2:end);
                
                Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = [IntPrb.W{1}.W(IntPrb.GridGamma,:),IntPrb.W{2}.W(IntPrb.GridGamma,:)]*cn;
                Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                 Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}.W(ExtPrb.GridGamma,:),ExtPrb.W{2}.W(ExtPrb.GridGamma,:)]*cn ;%- uinc;
                
                UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',ExtPrb.Scatterer.Eta);
                %UincParams  = struct('ScattererType','circle','r',PlrGrid.R);
                poUinc = Uinc(UincParams,ExtPrb.Scatterer.Phi,IncAng,k);
                
                
                %tmp = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
%                 ExtuInc = spalloc(PlrGrid.Nx,PlrGrid.Ny,length(ExtPrb.Scatterer.Nm));
%                 ExtuInc(ExtPrb.Scatterer.Outside) = tmp(ExtPrb.Scatterer.Outside);
                Extu = ExtPrb.P_Omega(Extxi,poUinc ); 
                
            else
                
              
                cn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*cn0ex )) ; 
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*cn0ex + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
                Extu = ExtPrb.P_Omega(Extxi,zeros(Nr,Nth-1));
                
                
                cn1 =( IntPrb.Q1 \ ( -IntPrb.Q0*cn0ex )) ;
                
                Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*cn0ex + IntPrb.W1(IntPrb.GridGamma,:)*cn1;
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
                ExtErr(n) =norm(tmp(:),inf);
                
                tmp = Intu(1:2:end,1:2:end)-Intu1(1:2:end,1:2:end);
                IntErr(n) =norm(tmp(:),inf);
                
                fprintf('b=%-7.2f,kex=%d,kin=%d,NBss0=%d,NBss1=%d,Nplr=%-5dx%d\t, Ncrt=%-5dx%d\t ExtErr=%d\t IntErr=%d\t time=%d\n',b, ...
				ExtWaveNumberParams.k ,IntWaveNumberParams.k ,Basis.NBss0,Basis.NBss1, Nr,Nth,Nx,Ny,full(ExtErr(n)),full(IntErr(n)),t);
            end
            
            Extu0=spalloc(Nr*2-1,Nth*2-2,nnz(Extu));
            Extu0(1:2:end,1:2:end)=Extu;
            
            Intu0=spalloc(Nx*2-1,Ny*2-1,nnz(Intu));
            Intu0(1:2:end,1:2:end)=Intu;

            if want2plot %|| n==nmax
                
                
                
                R = ones(size(PlrGrid.R)).*NaN;
                Th = ones(size(PlrGrid.R)).*NaN;
                Nm = ExtPrb.Scatterer.Nm;
                R(Nm)=PlrGrid.R(Nm); 
                Th(Nm)=PlrGrid.Theta(Nm);
                
                R(:,Nth)=R(:,Nth-1);
                Th(:,Nth)=2*pi;
                
                XExt = R .* cos(Th);
                YExt = R .* sin(Th);
                
                
               % UincParams  = struct('ScattererType','circle','r',PlrGrid.R);
               %tExtu=Extu; %+ Uinc(UincParams,PlrGrid.Theta,IncAng,k);
               %tExtu(ExtPrb.Scatterer.Mp)=0;%NaN;
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
                
                filename = sprintf('TFElpsAR%dRing%d%dkin%dkex%dgrd%d',fix(10*a/b),fix(10*r0),fix(10*r1),IntWaveNumberParams.k,ExtWaveNumberParams.k,Nr);
                ExScat = ExtPrb.Scatterer;
                IntScat = IntPrb.Scatterer;
                save(filename,'filename','a','b','XInt','XExt','YInt','YExt','tIntu','tExtu','Extu','Intu','Nm','Np','PlrGrid','CrtsGrid','ExScat','IntScat');

                continue
                
                th=0:0.0001:2*pi;
                
                figure(1)
                pcolor([XInt,XExt],[YInt,YExt],full(abs([tIntu,tExtu])))
                
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
                hold
                plot(a*cos(th),b*sin(th),'k','LineWidth',2)
                hold off
                
                figure(2)
                pcolor([XInt,XExt],[YInt,YExt],full(real([tIntu,tExtu])))
                
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
                hold
                plot(a*cos(th),b*sin(th),'k','LineWidth',2)
                hold off
                
                figure(3)
                pcolor([XInt,XExt],[YInt,YExt],full(imag([tIntu,tExtu])))
                
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
                 hold
                plot(a*cos(th),b*sin(th),'k','LineWidth',2)
                hold off
                
%                 tmp=abs([tIntu,tExtu]);
%                fprintf('min=%d \t max=%d', full(min(tmp(:))), full(max(tmp(:))));
%                 
              % saveas(figure(1),[filename 'abs.jpg'],'jpg')
              % saveas(figure(2),[filename 'real.jpg'],'jpg')
              % saveas(figure(3),[filename 'imag.jpg'],'jpg')
            end
            
            
        end
        
        ExtL=log2(ExtErr(1:end-1)./ExtErr(2:end))
        IntL=log2(IntErr(1:end-1)./IntErr(2:end))
    end
    end
end


function uinc = Uinc(Params,phi,IncAng,k)  
    
    if strcmpi(Params.ScattererType,'ellipse')
         x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
         y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        x = Params.r .* cos(phi);
        y = Params.r .* sin(phi);
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
    end
    
   
    
    uinc = Uinc(Params,phi,IncAng,k)  ;
    
    duinc = 1i .* k .*  uinc .* (dx.*cos(IncAng) + dy.*sin(IncAng));
 %   h = FocalDist*sqrt(sinh(eta).^2 + sin(phi).^2);
   % duinc = duinc./h;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



