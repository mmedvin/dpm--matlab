function RunSimpleTransReflAboutCircle
    
    want2plot = 0;
    
    r0 = 0.6;
    r1 = 2;
    
    x1=-0.8;xn=0.8;
    y1=-0.8;yn=0.8;
   
    R0 =0.7;
    NHR = 1.6;
    
    IncAng = 40;
    IncAng = IncAng*pi/180;    
    
    dk = 0;
    
%     Problem = 'Dirichlet'; % 'Dirichlet' or 'Neumann'   
%     fprintf('Solving %s defraction problem, comparing using grid convergance, data is PlaneWave, scatterer is circle\n',Problem);
     fprintf('RunSimpleTransReflAboutCircle \n');
    for k = 1 %[1, 5,20,25];%[1,3,5,10]%[1,5,10,15,20,25]
            
        UincParams  = struct('ScattererType','circle','r',R0);              
        f1      = @(phi) Uinc(UincParams,phi,IncAng,k+dk);
        dfdn    = @(phi) detaUinc(UincParams,phi,IncAng,k+dk);
            

        
    %    [cn0ex,cn1ex,M] = FourierCoeff(f1,dfdn);
    
        nmax=5;
        IntErr=zeros(1,nmax);
        ExtErr=zeros(1,nmax);
        
        WaveNumberHandle = @Tools.Coeffs.ConstantWaveNumber;
        ScattererHandle  = @Tools.Scatterer.PolarScatterer;
        ScattererParams  = struct('r0',R0,'ExpansionType',15, 'Stencil', 9);
        
        for n=0:nmax %run different grids
            tic
            %build grid
            
            p=3;%3;%1;
            Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
            Nx=2^(n+p)+1;	Ny=2^(n+p)+1;
           
            %BasisIndices        = -M:M;
             Basis = Tools.Basis.FourierBasis.BasisHelper(f1,dfdn,[1e-13,1e-13]);
            
            PlrGrid                = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
            			
            ExtPrb = Solvers.ExteriorSolver( struct(...
                     'Basis',Basis, ...
                     'Grid', PlrGrid, ...
                     'CoeffsHandle', WaveNumberHandle, ... 
                     'CoeffsParams', struct('k',k,'r0',NHR), ...
                     'ScattererHandle',ScattererHandle, ...
                     'ScattererParams', ScattererParams, ...
                     'CollectRhs',0, ... %i.e. no
                     'Extension', @Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension, ...
                     'ExtensionParams',[] ...
                     ));

            CrtsGrid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
            
            IntPrb = Solvers.InteriorHomoSolver( struct(...
                     'Basis',Basis, ...
                     'Grid', CrtsGrid, ...
                     'CoeffsHandle', WaveNumberHandle, ... 
                     'CoeffsParams', struct('k',k+dk,'r0',NHR), ...
                     'ScattererHandle',ScattererHandle, ...
                     'ScattererParams', ScattererParams, ...
                     'CollectRhs',1, ... %i.e. yes
                     'Extension', @Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension, ...
                     'ExtensionParams',[] , ...
                     'DiffOp'            , @Tools.DifferentialOps.HelmholtzOp, ...
                     'DiffOpParams'      , [] ...
                     ));
                 
            if 1
                
                UincParams  = struct('ScattererType','circle','r',ExtPrb.Scatterer.r);
                
                rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
                % rhs(numel(IntPrb.GridGamma)+1:end,1)= -Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
                uinc = Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
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
                
                UincParams  = struct('ScattererType','circle','r',PlrGrid.R);
                Uinc_ = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
                
                
                %tmp = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
%                 ExtuInc = spalloc(PlrGrid.Nx,PlrGrid.Ny,length(ExtPrb.Scatterer.Nm));
%                 ExtuInc(ExtPrb.Scatterer.Outside) = tmp(ExtPrb.Scatterer.Outside);
                Extu = ExtPrb.P_Omega(Extxi,Uinc_ ); 
                
            else
                
              
                cn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*cn0ex )) ; 
                
                Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*cn0ex + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
                Extu = ExtPrb.P_Omega(Extxi);
                
                
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
                
                fprintf('kex=%d,kin=%d,NBss0=%d, NBss1=%d,Nplr=%-5dx%d\t, Ncrt=%-5dx%d\t ExtErr=%d\t IntErr=%d\t time=%d\n',k,k+dk,Basis.NBss0,Basis.NBss1, Nr,Nth,Nx,Ny,full(ExtErr(n)),full(IntErr(n)),t);
            end
            
            Extu0=spalloc(Nr*2-1,Nth*2-2,nnz(Extu));
            Extu0(1:2:end,1:2:end)=Extu;
            
            Intu0=spalloc(Nx*2-1,Ny*2-1,nnz(Intu));
            Intu0(1:2:end,1:2:end)=Intu;

            if want2plot || n==nmax
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
                
%                 tmp=abs([tIntu,tExtu]);
%                fprintf('min=%d \t max=%d', full(min(tmp(:))), full(max(tmp(:))));
                
%                saveas(figure(1),'totfldkin15kout30abs.jpg','jpg')
%                saveas(figure(2),'totfldkin15kout30real.jpg','jpg')
%                saveas(figure(3),'totfldkin15kout30imag.jpg','jpg')
            end
            
            
        end
        
        ExtL=log2(ExtErr(1:end-1)./ExtErr(2:end))
        IntL=log2(IntErr(1:end-1)./IntErr(2:end))
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



