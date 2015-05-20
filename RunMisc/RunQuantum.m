function RunQuantum
    
    want2plot = 0;
 
    %ellipse params
    a=1;     b=a/2;
    FocalDistance = sqrt(a^2-b^2); assert(isreal(FocalDistance ));
    Eta0 = acosh(a/FocalDistance);
    
    %exterior problem ring params
    r0 = 0.3;   assert(r0<b);
    r1 = 2.2;   assert(r1>a);
    
    %interior problem rectangular params
    x1=-1.2; xn=1.2;
    y1=-0.7; yn=0.7;
   
    %incident angle
    IncAng = 40*pi/180;    

    Vin = 2;
    Vex = 1;  

%------------------------------------------------------------------------------------------------------------------------------------        
   IntOp=0; %non zero value - use helmholtz homogeneous solver     
%-----------------------TMP-------------------------------------------
    ExtVAddParams = struct('k',Vex); 
    IntVAddParams = struct('k',Vin,'r0',1.6); 
%-----------------------TMP-------------------------------------------
        
    UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',Eta0, 'Vark',true); % ,false);%
    f1      = @(phi) Uinc(UincParams,phi,IncAng,Vin);
    dfdn    = @(phi) detaUinc(UincParams,phi,IncAng,Vex);
        
    % if it doesn't work well automatically - ask me about the other two options
    Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn); % estimate how many basis functions to use automatically
    %Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn, 10^-12); % estimate how many basis functions to use given error 
    %Basis =Tools.Basis.FourierBasis.BasisHelper(f1,dfdn, 50); %  use given number of basis functions to use given error 
    
    fprintf('Quantum Prb: ELPS: M=%d;  a=%-4.2f, b=%-4.2f; Ring: r0=%-4.2f, r1=%-4.2f; Rec:(x1,xn)x((y1,yn)=(%-4.2f,%-4.2f)x(%-4.2f,%-4.2f); IncAng=%-4.2f^deg; Vin=%-4.2f, Vex=%-4.2f \n',...
            Basis.M, a,b,r0,r1,x1,xn,y1,yn,IncAng*180/pi,Vin, Vex);

    
    
    IntErr=0;ExtErr=0;
    LastIntErr=0;LastExtErr=0;
    
    %run different grids to meausure the grid convergence
    nmax=3;
    for n=1:nmax 
        tic
        
        %build grid        
        p=4;%6;
        Nr=2^(n+p)+1;	Nth=2^(n+p)+1;
        Nx=2^(n+p)+1;	Ny=2^(n+p)+1;
       
        %---------------------------------------------------------------------        
        ScattererHandle  = @Tools.Scatterer.EllipticScatterer;        
        %-----------------------EXTERIOR--------------------------------------        
        ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',25, 'Stencil', 9);
        PlrGrid                = Tools.Grid.PolarGrids(r0,r1,Nr,Nth);
        WaveNumberClsHandle = @Tools.Coeffs.ConstantWaveNumber;
        CollectRhs=0;
        
        ExtPrb = Solvers.ExteriorSolver ...
            (Basis,PlrGrid,WaveNumberClsHandle,ExtVAddParams,ScattererHandle,ScattererParams,CollectRhs);
       
        %-----------------------INTERIOR--------------------------------------
        CrtsGrid                = Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny);
        
        WaveNumberClsHandle = @Tools.Coeffs.WaveNumberElliptical;
        
        if IntOp
            ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',25, 'Stencil', 9);
            CollectRhs=1;
            IntPrb = Solvers.InteriorHomoSolver ...
                (Basis,CrtsGrid,WaveNumberClsHandle,IntVAddParams,ScattererHandle,ScattererParams,CollectRhs);
        else
            ScattererParams  = struct('Eta0',Eta0,'FocalDistance',FocalDistance,'ExpansionType',35, 'Stencil', 13);
            LinearSolverType = 0;
            DiffOp = @Tools.DifferentialOps.LapOp4OrdrVarCoeffBCinRhs;
            
            InteriorCoeffsHandle = @Tools.Coeffs.LaplaceCoeffsEllps2;
            InteriorCoeffsParams = struct('ca',1,'da',0,'ea',0,'WithB',1,'cb',1,'db',1,'eb',-2,'sigma',0);
            %InteriorCoeffsParams = struct('ca',111,'da',0,'ea',0,'WithB',1,'cb',55,'db',0,'eb',-2,'sigma',0);
            
            DiffOpParamsInt = struct('BC_x1',0,'BC_xn', 0,'BC_y1',0,'BC_yn',0, 'LinearSolverType', LinearSolverType, 'Order',4);
            
            
            IntPrb = Solvers.InteriorHomoLaplacianSolver ...
                (Basis,CrtsGrid,InteriorCoeffsHandle,InteriorCoeffsParams,ScattererHandle,ScattererParams,CollectRhs,DiffOp,DiffOpParamsInt);
            
        end
        if 1 % solve it
            UincParams  = struct('ScattererType','ellipse','FocalDistance',FocalDistance,'eta',ExtPrb.Scatterer.eta, 'Vark',false);
            
            rhs = zeros(numel(ExtPrb.GridGamma) + numel(IntPrb.GridGamma),1);
            uinc = Uinc(UincParams,ExtPrb.Scatterer.phi,IncAng,ExtVAddParams.k);
            rhs(numel(IntPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
            
            cn = [ IntPrb.Q; ExtPrb.Q ] \ rhs;
            
            Intxi = spalloc(Nx,Ny   ,length(IntPrb.GridGamma));
            Intxi(IntPrb.GridGamma) = IntPrb.W(IntPrb.GridGamma,:)*cn;
            Intu = IntPrb.P_Omega(Intxi);
            
            Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
            Extxi(ExtPrb.GridGamma) = ExtPrb.W(ExtPrb.GridGamma,:)*cn ;
            
            UincParams  = struct('ScattererType','circle','r',PlrGrid.R, 'Vark',false);
            Uinc = Uinc(UincParams,PlrGrid.Theta,IncAng,ExtVAddParams.k);
            
            Extu = ExtPrb.P_Omega(Extxi,Uinc );            
            
        else  % for debug only
            cn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*cn0ex )) ;
            
            Extxi = spalloc(Nr,Nth-1,length(ExtPrb.GridGamma));
            Extxi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*cn0ex + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
            Extu = ExtPrb.P_Omega(Extxi,zeros(Nr,Nth-1));
            
            cn1 =( IntPrb.Q1 \ ( -IntPrb.Q0*cn0ex )) ;
            
            Intxi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
            Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*cn0ex + IntPrb.W1(IntPrb.GridGamma,:)*cn1;
            Intu = IntPrb.P_Omega(Intxi);                       
        end
        
        t=toc;
        
%-----------------------------------------------------------------------------
% Comparison
%-----------------------------------------------------------------------------
        if n > 1 %nothing to compare with yet
            Extu1= spalloc(Nr,Nth-1,nnz(Extu0));
            Extu1(ExtPrb.Scatterer.Nm) = Extu0(ExtPrb.Scatterer.Nm);
            
            Intu1= spalloc(Nx,Ny,nnz(Intu0));
            Intu1(IntPrb.Np) = Intu0(IntPrb.Np);
            
            tmp = Extu(1:2:end,1:2:end)-Extu1(1:2:end,1:2:end);
            ExtErr =norm(tmp(:),inf);
            
            tmp = Intu(1:2:end,1:2:end)-Intu1(1:2:end,1:2:end);
            IntErr =norm(tmp(:),inf);
            
            fprintf('Nplr=%-5dx%-5d\t, Ncrt=%-5dx%-5d\t Err(Int|Ext)=%-7.6d|%-7.6d\t Rate(Int|Ext)=%-5.4f|%-5.4f\t time=%d\n', Nr,Nth,Nx,Ny,...
                IntErr,   ExtErr,   log2(LastIntErr/IntErr),     log2(LastExtErr/ExtErr),t);
            
            LastIntErr = IntErr;
            LastExtErr = ExtErr;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



