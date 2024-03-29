function RunTransReflAboutElps
    
    want2plot = 0;

    IntErr=[];ExtErr=[];
    
    a=1;%2.5;
    b=0.8;
    
    r0 = 0.7*b;
    r1 = 1.8*a;
    
    x1=-1.2;xn=1.2;
    y1=-1;  yn=1;
   
  %  R0 =0.7;
   
    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
  
    IncAng = 40;
    IncAng = IncAng*pi/180;    
        
%     Problem = 'Dirichlet'; % 'Dirichlet' or 'Neumann'   
%     fprintf('Solving %s defraction problem, comparing using grid convergance, data is PlaneWave, scatterer is circle\n',Problem);
     fprintf('RunTransReflAboutElps \n');
kin = [15];
kex = [1];

    for ki = 1 %[1, 5,20,25];%[1,3,5,10]%[1,5,10,15,20,25]
    ExtWaveNumberParams = struct('k',kex(ki),'r0',1.6);
    IntWaveNumberParams = struct('k',kin(ki),'r0',1.6);
            
        UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',Eta0, 'Vark',true); % ,false);%
        f1      = @(phi) full(Uinc(UincParams,phi,IncAng,kin(ki)));
        dfdn    = @(phi) full(detaUinc(UincParams,phi,IncAng,kin(ki)));
            

        Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn);
    
        nmax=4;
        for n=1:nmax %run different grids
            tic
            %build grid
            
            p=5;%6;%3;%1;
            Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
            Nx=2^(n+p)+1;	Ny=2^(n+p)+1;
           
%             BasisIndices        = -M:M;
            PlrGrid                = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
                      
            ScattererHandle = @Tools.Scatterer.EllipticScatterer;
            ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'Stencil',9);
            			
            ExtPrb = Solvers.ExteriorSolver( struct(...
                     'Basis',Basis, ...
                     'Grid', PlrGrid, ...
                     'CoeffsHandle', @Tools.Coeffs.ConstantWaveNumber, ... 
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
                     'CoeffsHandle', @Tools.Coeffs.WaveNumberElliptical, ... 
                     'CoeffsParams', IntWaveNumberParams, ...
                     'ScattererHandle',ScattererHandle, ...
                     'ScattererParams', ScattererParams, ...
                     'CollectRhs',1, ... %i.e. yes
                     'Extension', @Tools.Extensions.TwoTupleExtension, ...
                     'ExtensionParams',[] , ...
                     'DiffOp'            , @Tools.DifferentialOps.HelmholtzOp, ...
                     'DiffOpParams'      , [] ...
                     ));
         
            if 1
         
                UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',ExtPrb.Scatterer.eta, 'Vark',false);
                
                rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
                % rhs(numel(IntPrb.GridGamma)+1:end,1)= -Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
                uinc = Uinc(UincParams,ExtPrb.Scatterer.phi,IncAng,ExtWaveNumberParams.k);
                rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
                
                cn = [ IntPrb.Q{1},IntPrb.Q{2} ; ExtPrb.Q{1},ExtPrb.Q{2} ] \ rhs;
                
                %             cn0 = cn(1:2*M+1);
                %             cn1 = cn(2*M+2:end);
                
                Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = [IntPrb.W{1}(IntPrb.GridGamma,:),IntPrb.W{2}(IntPrb.GridGamma,:)]*cn;
                Intu = IntPrb.P_Omega(Intxi);
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}(ExtPrb.GridGamma,:),ExtPrb.W{2}(ExtPrb.GridGamma,:)]*cn ;%- uinc;
                %Extxi(ExtPrb.GridGamma) = [ExtPrb.W(ExtPrb.GridGamma,:),uinc]*cn;
                
                %UincParams  =
                %struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',PlrGrid.R);???
                UincParams  = struct('ScattererType','circle','r',PlrGrid.R, 'Vark',false);
                Uinc_ = Uinc(UincParams,PlrGrid.Theta,IncAng,ExtWaveNumberParams.k);
                
                
                %tmp = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
%                 ExtuInc = spalloc(PlrGrid.Nx,PlrGrid.Ny,length(ExtPrb.Scatterer.Nm));
%                 ExtuInc(ExtPrb.Scatterer.Outside) = tmp(ExtPrb.Scatterer.Outside);
                Extu = ExtPrb.P_Omega(Extxi,Uinc_ ); 
                
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
                
                fprintf('kex=%d,kin=%d,NBss0=%d,NBss1=%d,Nplr=%-5dx%d\t, Ncrt=%-5dx%d\t ExtErr=%d\t IntErr=%d\t time=%d\n',...
                    ExtWaveNumberParams.k,IntWaveNumberParams.k,Basis.NBss0,Basis.NBss1, Nr,Nth,Nx,Ny,full(ExtErr(n)),full(IntErr(n)),t);
            end
            
            Extu0=spalloc(Nr*2-1,Nth*2-2,nnz(Extu));
            Extu0(1:2:end,1:2:end)=Extu;
            
            Intu0=spalloc(Nx*2-1,Ny*2-1,nnz(Intu));
            Intu0(1:2:end,1:2:end)=Intu;

            if want2plot %|| n==nmax
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
                close all
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
                close all
                
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
                close all
                
                tmp=abs([tIntu,tExtu]);
               fprintf('min=%d \t max=%d', full(min(tmp(:))), full(max(tmp(:))));
                
               saveas(figure(1),'totfldVarKkin50kout15abs.jpg','jpg')
               saveas(figure(2),'totfldVarKkin50kout15real.jpg','jpg')
               saveas(figure(3),'totfldVarKkin50kout15imag.jpg','jpg')
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
    end
 
    if Params.Vark
        IntWaveNumberAddParams = struct('k',k0,'r0',1.6);
        Scat = struct('FocalDistance',Params.FocalDistance,'Eta',Params.eta,'Phi',phi);
        WN = Tools.Coeffs.WaveNumberElliptical(Scat, IntWaveNumberAddParams );
        k=WN.k;
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



function [cn0,cn1,M] = FourierCoeff(f,dfdr)
cn0 =   quadgk(f,-pi,pi);
cn1  =   quadgk(dfdr,-pi,pi);
j=0;
notstop=true;
err=2*pi*(1e-10);
while notstop
    j=j+1;
    g = @(x) f(x).*exp(-1i*j*x);
    tcn0p =   quadgk(g,-pi,pi);
    
    g = @(x) f(x).*exp(-1i*(-j)*x);
    tcn0m =   quadgk(g,-pi,pi);
    
    dg = @(x) dfdr(x).*exp(-1i*j*x);
    tcn1p  =   quadgk(dg,-pi,pi);
    
    dg = @(x) dfdr(x).*exp(-1i*(-j)*x);
    tcn1m  =   quadgk(dg,-pi,pi);
    
    
    cn0 = [tcn0m; cn0 ;tcn0p];
    cn1 = [tcn1m; cn1 ;tcn1p];
    
   % notstop  = abs(cn0(1))>err;
    notstop  = abs(cn0(1))>err | abs(cn0(2))>err | ...
                    abs(cn0(end))>err | abs(cn0(end-1))>err;
   
end
cn0=cn0/(2*pi);
cn1=cn1/(2*pi);
M=j;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



