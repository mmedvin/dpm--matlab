function RunTransReflAboutSubmarine
 
LastKnownVersionOfTools = '2.0.0.0';
try
    ver = Tools.Version();
    if ~strcmp(ver,LastKnownVersionOfTools)
        error('MDP:wrong version of Tools, expected version %s, found version %s',LastKnownVersionOfTools,ver);
    end
catch err
    if strcmp(err.identifier, 'MATLAB:undefinedVarOrClass')
path
        error('MDP: please add parent folder to the path');
    else%if strcmp(err.identifier,'MDP:wrong version of Tools')
        sprintf(err.message);
        rethrow(err);
    end
end   
    want2plot = 0;
    AR=2;
   a=1.8;    b=a/AR;
   r0=0.3; r1=3;
    
    x1=-3;xn=3;
    y1=-3;  yn=3;
   
    ellipse = struct('a',a,'b',b);
    tower = struct('c',4,'p',20);
%	tower = struct('c',0,'p',10);
    
%     FocalDistance = sqrt(a^2-b^2);
%     Eta0 = acosh(a/FocalDistance);
  
    IncAng = 40;
    IncAng = IncAng*pi/180;    
    
%     dk = 30;
    
    ScatType = 'submarine'; %'ellipse' or 'circle' or 'submarine'
%     BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
    
%fprintf('Trans/Refl problem about Submarine of FD=%d, ,Eta0=%d, a=%d, b=%d,c=%d,p=%d, AR=%d ,r0=%d,r1=%d,x=+/-%d,y+/-=%d\n',FocalDistance , Eta0 , ellipse.a,ellipse.b,tower.c,tower.p,a/b,r0,r1,xn,yn);
   fprintf('Trans/Refl problem about Submarine of FD=%d, ,Eta0=%d, a=%d, b=%d,c=%d,p=%d, AR=%d ,r0=%d,r1=%d,x=+/-%d,y+/-=%d\n',sqrt(a^2-b^2) , acosh(a/sqrt(a^2-b^2)) , ellipse.a,ellipse.b,tower.c,tower.p,ellipse.a/ellipse.b,r0,r1,xn,yn);
     

kex = [1,5,10];
kin = kex*3;

    for ki =1% [1, 5,20,25];%[1,3,5,10]%[1,5,10,15,20,25]
    ExtWaveNumberAddParams = kex(ki);
    IntWaveNumberAddParams = kin(ki);%struct('k0',kin(ki),'r0',1.6);
     k=ExtWaveNumberAddParams;
    
        
      if strcmpi(ScatType,'ellipse')
        UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',Eta0, 'Vark',true); % ,false);%
    elseif strcmpi(ScatType,'circle')
      %  ExParams  = struct('ScattererType','circle','r',R0, 'HankelIndex', HankelIndex,'HankelType',HankelType);
    elseif strcmpi(ScatType,'submarine')       
        UincParams  = struct('ScattererType','submarine','ellipse',ellipse,'tower',tower, 'Vark',true);
    end

    
        
        f1      = @(phi) Uinc(UincParams,phi,IncAng,kin(ki));
        dfdn    = @(phi) detaUinc(UincParams,phi,IncAng,kin(ki));
            

        Basis =Tools.Basis.FourierBasis.BasisHelper(f1,@sin);%dfdn);
        %[cn0ex,cn1ex,M] = FourierCoeff(f1,dfdn);
    
        nmax=5;%8
        for n=1:nmax %run different grids
            tic
            %build grid
            
            p=4;%3;%1;
            Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
            Nx=2^(n+p)+1;	Ny=2^(n+p)+1;
           
%             BasisIndices        = -M:M;
            PlrGrid                = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
            WaveNumberClsHandle = @Tools.WaveNumber.ConstantWaveNumber;
                      
            ScattererClsHandle  = @Tools.Scatterer.SubmarineScatterer;
            ScattererAddParams  = struct('ellipse',ellipse,'tower',tower,'ExpansionType',25);
            
            ExtPrb = Solvers.ExteriorSolver ...
                (Basis,PlrGrid,WaveNumberClsHandle,ExtWaveNumberAddParams,ScattererClsHandle,ScattererAddParams);

            CrtsGrid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);

            WaveNumberClsHandle = @Tools.WaveNumber.ConstantWaveNumber;%WaveNumberElliptical;
            

            %IntWaveNumberAddParams = k+dk;           
            
            IntPrb = Solvers.InteriorHomoSolver ...
                (Basis,CrtsGrid,WaveNumberClsHandle,IntWaveNumberAddParams,ScattererClsHandle,ScattererAddParams);
         
            if 0
                IntQ = IntPrb.Q;
                ExtQ = ExtPrb.Q;%[ExtPrb.Q0,-ExtPrb.Q1]; %
                
                %%%%%%%%%%%%%%%%%5
                if strcmpi(ScatType,'ellipse')
                    UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',ExtPrb.Scatterer.eta, 'Vark',false);                    
                elseif strcmpi(ScatType,'circle')
%                     ExParams3 =ExParams;
%                     ExParams3.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
%                     exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
                elseif strcmpi(ScatType,'submarine')
                    UincParams  = struct('ScattererType','submarine','ellipse',ellipse,'tower',tower, 'Vark',true);                    
                    UincParams.r = ExtPrb.Scatterer.r;
                   UincParams.Vark=false;
                end
                uinc = Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
                %%%%%%%%%%%%%%%%%%%%%%
                             
                rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
                rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
                
                cn = [ IntQ ; ExtQ ] \ rhs;
                
                %             cn0 = cn(1:2*M+1);
                %             cn1 = cn(2*M+2:end);
                
                Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*cn;
                Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = ExtPrb.W(ExtPrb.GridGamma,:)*cn ;%- uinc;
                %Extxi(ExtPrb.GridGamma) = [ExtPrb.W(ExtPrb.GridGamma,:),uinc]*cn;
                
                %UincParams  =
                %struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',PlrGrid.R);???
               % UincParams  = struct('ScattererType','circle','r',PlrGrid.R, 'Vark',false);
               UincParams  = struct('ScattererType','submarine','ellipse',ellipse,'tower',tower, 'Vark',false); 
               UincParams.r=PlrGrid.R;
               Uinc = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
                
                
                %tmp = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
%                 ExtuInc = spalloc(PlrGrid.Nx,PlrGrid.Ny,length(ExtPrb.Scatterer.Nm));
%                 ExtuInc(ExtPrb.Scatterer.Outside) = tmp(ExtPrb.Scatterer.Outside);
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
                ExtErr(n) =norm(tmp(:),inf);
                
                tmp = Intu(1:2:end,1:2:end)-Intu1(1:2:end,1:2:end);
                IntErr(n) =norm(tmp(:),inf);
                
                fprintf('kex=%d,kin=%d,M=%d,Nplr=%-5dx%d\t, Ncrt=%-5dx%d\t ExtErr=%d\t IntErr=%d\t time=%d\n',kex(ki),kin(ki),Basis.M, Nr,Nth,Nx,Ny,full(ExtErr(n)),full(IntErr(n)),t);
            end
            
            Extu0=spalloc(Nr*2-1,Nth*2-2,nnz(Extu));
            Extu0(1:2:end,1:2:end)=Extu;
            
            Intu0=spalloc(Nx*2-1,Ny*2-1,nnz(Intu));
            Intu0(1:2:end,1:2:end)=Intu;

            if want2plot && n==nmax
                figure
                XExt = PlrGrid.R .* cos(PlrGrid.Theta);
                YExt = PlrGrid.R .* sin(PlrGrid.Theta);
                XExt(ExtPrb.Scatterer.Mp)=NaN;
                YExt(ExtPrb.Scatterer.Mp)=NaN;
                
                
               % UincParams  = struct('ScattererType','circle','r',PlrGrid.R);
                tExtu=Extu; %+ Uinc(UincParams,PlrGrid.Theta,IncAng,k);
                tExtu(ExtPrb.Scatterer.Mp)=NaN;
                
                tIntu = Intu;
                tIntu(IntPrb.Scatterer.Mm)=NaN;
                
                XInt = CrtsGrid.X;
                YInt = CrtsGrid.Y;
                XInt(IntPrb.Scatterer.Mm)=NaN;
                YInt(IntPrb.Scatterer.Mm)=NaN;
                
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
                saveas(figure(1),'totfldElps_a1.6b0.8r00.6r12ang70VarKkin45kout15abs.jpg','jpg')
               % close all
                figure
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
                saveas(figure(1),'totfldElps_a1.6b0.8r00.6r12ang70VarKkin45kout15real.jpg','jpg')
                %close all
                
                figure
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
                saveas(figure(1),'totfldElps_a1.6b0.8r00.6r12ang70VarKkin45kout15imag.jpg','jpg')
                %close all
                
                tmp=abs([tIntu,tExtu]);
               fprintf('min=%d \t max=%d', full(min(tmp(:))), full(max(tmp(:))));
                
%                saveas(figure(1),'totfldVarKkin50kout15abs.jpg','jpg')
%                saveas(figure(2),'totfldVarKkin50kout15real.jpg','jpg')
%                saveas(figure(3),'totfldVarKkin50kout15imag.jpg','jpg')
            end
            
            
        end
        
        ExtL=log2(ExtErr(1:end-1)./ExtErr(2:end))
        IntL=log2(IntErr(1:end-1)./IntErr(2:end))
    end
    
end


function uinc = Uinc(Params,phi,IncAng,k0)  
    
    if strcmpi(Params.ScattererType,'ellipse')
         x = Params.FocalDistance * cosh(Params.eta) .* cos(phi);
         y = Params.FocalDistance * sinh(Params.eta) .* sin(phi);
    elseif strcmpi(Params.ScattererType,'circle')
        x = Params.r .* cos(phi);
        y = Params.r .* sin(phi);
    elseif strcmpi(Params.ScattererType,'submarine')        
        try
            r=Params.r  ;
        catch
            e = (cos(phi).^2/Params.ellipse.a^2 + sin(phi).^2/Params.ellipse.b^2);
            r = sqrt((1 + Params.tower.c*sin(phi).^Params.tower.p) ./ e);
        end
                
        x = r .* cos(phi);
        y = r .* sin(phi);
    end
 
    if Params.Vark        
        k=k0;%Tools.WaveNumber.ConstantWaveNumber. StarShaped.kkn(ellipse,tower,theta,k0,r0);
    else
        k=k0;
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



