function RunTransReflAboutStarShapedBody
    
    want2plot = 0;

    IntErr=[];ExtErr=[];
    
    a=1;%2.5;
    b=a/2;%0.8;
    
    %r0 = 0.7*b;
    %r1 = 1.8*a;
    
r0=0.3;
r1=2.2;

    x1=-2;xn=2;
    y1=-2;  yn=2;
   
  %  R0 =0.7;
   
    FocalDistance = sqrt(a^2-b^2);
    Eta0 = acosh(a/FocalDistance);
  
    IncAng = 40;
    IncAng = IncAng*pi/180;    
    
    dk = 30;
    
    %doesn't expected to work Parameterization  = Tools.Parameterizations.ParametricHeart(struct('a',13/16,'b',-5/16,'c',-2/16,'d',-1/16,'e',1,'p',3));
    Parameterization  = Tools.Parameterizations.ParametricEllipse(struct('a',a,'b',b));
    %Parameterization  = Tools.Parameterizations.ParametricKite(struct('a',1,'b',.65*2,'c',1.5));
    %Parameterization  = Tools.Parameterizations.ParametricSubmarine(struct('a',1,'b',1/2,'c',0,'p',200));
    %Parameterization  = Tools.Parameterizations.ParametricStar();
    
    
%     Problem = 'Dirichlet'; % 'Dirichlet' or 'Neumann'   
%     fprintf('Solving %s defraction problem, comparing using grid convergance, data is PlaneWave, scatterer is circle\n',Problem);
     
kin = [3,15,30];
kex = [1 ,5, 10];

    for ki = 1:3
        
        ErrIntPre = 0; 
        ErrExtPre = 0;
        
    ExtWaveNumberAddParams = struct('k',kex(ki),'r0',1.6);
    IntWaveNumberAddParams = struct('k',kin(ki),'r0',1.6);
    k=ExtWaveNumberAddParams.k;
    
        
        %UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',Eta0, 'Vark',true); % ,false);%
        UincParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization);
        f1      = @(phi) Uinc(UincParams,phi,IncAng,kin(ki));
        dfdn    = @(phi) detaUinc(UincParams,phi,IncAng,kin(ki));
            

        Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn);
        %[cn0ex,cn1ex,M] = FourierCoeff(f1,dfdn);
    
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
                      
            %IntScattererClsHandle  = @Tools.Scatterer.EllipticScatterer;%Interior
            %ExtScattererClsHandle  = @Tools.Scatterer.EllipticScatterer;%Exterior
            %ScattererAddParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance);
            
            ScattererHandle  = @Tools.Scatterer.StarShapedScatterer;
            ScattererParams  = UincParams;
            
			CollectRhs=0;
			
            ExtPrb = Solvers.ExteriorSolver ...
                (Basis,PlrGrid,WaveNumberHandle,ExtWaveNumberAddParams,ScattererHandle,ScattererParams,CollectRhs);

            CrtsGrid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);

            %WaveNumberHandle = @Tools.Coeffs.WaveNumberElliptical;
            

            %IntWaveNumberAddParams = k+dk;           
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
                uinc = Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
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
                %UincParams  = struct('ScattererType','circle','r',PlrGrid.R, 'Vark',false);
                UincParams = struct('ScattererType','StarShapedScatterer','Parameterization',Parameterization,'r',PlrGrid.R);
                Uinc = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
                
                
                %tmp = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
%                 ExtuInc = spalloc(PlrGrid.Nx,PlrGrid.Ny,length(ExtPrb.Scatterer.Nm));
%                 ExtuInc(ExtPrb.Scatterer.Outside) = tmp(ExtPrb.Scatterer.Outside);
                Extu = ExtPrb.P_Omega(Extxi,Uinc ); 
                
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
                ErrExt =norm(tmp(:),inf);
                
                tmp = Intu(1:2:end,1:2:end)-Intu1(1:2:end,1:2:end);
                ErrInt =norm(tmp(:),inf);
                
               % fprintf('kex=%d,kin=%d,M=%d,Nplr=%-5dx%d\t, Ncrt=%-5dx%d\t ExtErr=%d\t IntErr=%d\t time=%d\n',k,k+dk,Basis.M, Nr,Nth,Nx,Ny,full(ExtErr(n)),full(IntErr(n)),t);
               fprintf('kex=%d,kin=%d,M=%d,N=%-10dx%-10d\t ErrExt=%d\t rate=%-6.2f  ErrInt=%d\t rate=%-6.2f\t time=%d\n',...
			kex(ki),kin(ki),Basis.M, Nr,Nth,ErrExt,log2(ErrExtPre/ErrExt),ErrInt,log2(ErrIntPre/ErrInt),t);

       		ErrIntPre = ErrInt ; 
        	ErrExtPre = ErrExt ;

            end
            
            Extu0=spalloc(Nr*2-1,Nth*2-2,nnz(Extu));
            Extu0(1:2:end,1:2:end)=Extu;
            
            Intu0=spalloc(Nx*2-1,Ny*2-1,nnz(Intu));
            Intu0(1:2:end,1:2:end)=Intu;

            if 0 %want2plot %|| n==nmax
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
	   %duinc = duinc./h;
   end
   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



