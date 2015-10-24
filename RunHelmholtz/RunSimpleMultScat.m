function RunSimpleMultScat
    nmax=5;
    want2plot = 0;
    
    Boundaries=struct( ...
                            'Interior', struct('x1',-0.7, 'xn', 0.7, 'y1',-0.7,'yn',0.7), ...
                            'Medium'  , struct('x1',-1.2, 'xn', 1.2, 'y1',-1.2,'yn',1.2), ...
                            'Exterior', struct('r1', 0.8, 'rn', 1.5) ...
                      ); 
        
    R0 = 0.5;
    R1 = 1;
    
    NHR = struct('Interior', 1.6, 'Medium', 1.6);
    
    IncAng = Tools.Common.Angles(40,'degrees');
        
    ScatType = 'circle';
    BType    = 'Fourier';
    
    K = { struct( 'Interior',10,'Medium', 5, 'Exterior',1) ,  struct( 'Interior',1,'Medium',2, 'Exterior',1), struct( 'Interior',5,'Medium', 5, 'Exterior',5)};

    fprintf('RunSimpleMultScat, Incident Angle = %d deg \n', IncAng.Degrees);

    for k = 1:3 %[1, 5,20,25];%[1,3,5,10]%[1,5,10,15,20,25]
            
        if strcmpi(ScatType,'circle')
            UincParams{1}  = struct('ScattererType','circle', 'r', R0);
            UincParams{2}  = struct('ScattererType','circle', 'r', R1);            
        end
        
        f      = @(phi,n) Uinc    (UincParams{n}, phi, IncAng.Radians, K{k}.Interior);
        dfdn   = @(phi,n) detaUinc(UincParams{n}, phi, IncAng.Radians, K{k}.Interior);
            
        g      = @(phi,n) Uinc    (UincParams{n}, phi, IncAng.Radians, K{k}.Medium);
        dgdn   = @(phi,n) detaUinc(UincParams{n}, phi, IncAng.Radians, K{k}.Medium);

        h      = @(phi,n) Uinc    (UincParams{n}, phi, IncAng.Radians, K{k}.Exterior);
        dhdn   = @(phi,n) detaUinc(UincParams{n}, phi, IncAng.Radians, K{k}.Exterior);
        
        if strcmpi(BType,'Fourier')
            %Basis = struct( 'Interior', Tools.Basis.FourierBasis.BasisHelper(@(th) f(th,1) - g(th,1),@(th) dfdn(th,1) - dgdn(th,1)), ...
            %                'Exterior', Tools.Basis.FourierBasis.BasisHelper(@(th) g(th,1) - h(th,1),@(th) dgdn(th,1) - dhdn(th,1)) ...
            %    );
            Basis = struct( 'Interior', Tools.Basis.FourierBasis.BasisHelper(@(th) g(th,1) ,@(th) dgdn(th,1)), ...
                            'Exterior', Tools.Basis.FourierBasis.BasisHelper(@(th) g(th,2) ,@(th) dgdn(th,2)) ...
                );
            
        end
        
        TotErrPre=0;IntErrPre=0; MedErrPre=0; ExtErrPre=0;
        
        %WaveNumberHandle = struct('Interior', @Tools.Coeffs.ConstantWaveNumber               , 'Medium', @Tools.Coeffs.WaveNumberPolarR          , 'Exterior', @Tools.Coeffs.ConstantWaveNumber  );
        %WaveNumberHandle = struct('Interior', @Tools.Coeffs.WaveNumberPolarR               , 'Medium', @Tools.Coeffs.ConstantWaveNumber          , 'Exterior', @Tools.Coeffs.ConstantWaveNumber  );
        WaveNumberHandle = struct('Interior', @Tools.Coeffs.WaveNumberPolarR               , 'Medium', @Tools.Coeffs.WaveNumberPolarR          , 'Exterior', @Tools.Coeffs.ConstantWaveNumber  );
        %WaveNumberHandle = struct('Interior', @Tools.Coeffs.ConstantWaveNumber              , 'Medium', @Tools.Coeffs.ConstantWaveNumber       , 'Exterior', @Tools.Coeffs.ConstantWaveNumber  );
        WaveNumberParams = struct('Interior', struct('k',K{k}.Interior,'r0',NHR.Interior)  , 'Medium', struct('k',K{k}.Medium,'r0',NHR.Medium) , 'Exterior', struct('k',K{k}.Exterior)         );

        ScattererHandle  = struct('Interior', @Tools.Scatterer.PolarScatterer              , 'Medium', @Tools.Scatterer.NestedPolarScatterer   , 'Exterior', @Tools.Scatterer.PolarScatterer  );
        ScattererParams  = struct('Interior', struct('r0', R0, 'Stencil', 9)               , 'Medium', struct('r0',R0,'r1',R1, 'Stencil', 9)   , 'Exterior', struct('r0', R1, 'Stencil', 9)   );
                
        Extension       = @Tools.Extensions.NestedExtension;        
        ExtentionParams = struct('IntExtension',@Tools.Extensions.EBPolarHelmholtz5OrderExtension,'ExtExtension',@Tools.Extensions.EBPolarHelmholtz5OrderExtension);
        
        fprintf('Int xy: (%-2.2f,%-2.2f)x(%-2.2f,%-2.2f), Med xy: (%-2.2f,%-2.2f)x(%-2.2f,%-2.2f), Ext (%-2.2f,%-2.2f)\n',  Boundaries.Interior.x1, Boundaries.Interior.xn, Boundaries.Interior.y1, Boundaries.Interior.yn, ...
                                                                                                 Boundaries.Medium.x1  , Boundaries.Medium.xn  , Boundaries.Medium.y1  , Boundaries.Medium.yn  , ...
                                                                                                 Boundaries.Exterior.r1,Boundaries.Exterior.rn);
        fprintf('NBss: Int0=%d, Int1=%d, Ext0=%d, Ext1=%d\n',Basis.Interior.NBss0,Basis.Interior.NBss1,Basis.Exterior.NBss0,Basis.Exterior.NBss1);
        fprintf('k: Int: %s,%d, Med: %s,%d, Ext %s, %d \n', func2str(WaveNumberHandle.Interior), K{k}.Interior, ...
                                                            func2str(WaveNumberHandle.Medium)  , K{k}.Medium, ...
                                                            func2str(WaveNumberHandle.Exterior), K{k}.Exterior);
        fprintf('Scat: %s, %s, %s, R0=%-2.2f, R1=%-2.2f\n', func2str(ScattererHandle.Interior ), func2str(ScattererHandle.Medium) , func2str(ScattererHandle.Exterior ) ,R0,R1      );
        
        for n=0:nmax %run different grids
            tic
            %build grid
            
            p=4;%3;%1;
            N = struct('Interior',struct('Nx',2^(n+p)+1,'Ny',2^(n+p)+1),'Medium',struct('Nx',2^(n+p)+1,'Ny',2^(n+p)+1),'Exterior',struct('Nr', 2^(n+p)+1,'Nth', 2^(n+p)+1));           
                                     
            Grids.Exterior = Tools.Grid.PolarGrids(Boundaries.Exterior.r1,Boundaries.Exterior.rn,N.Exterior.Nr,N.Exterior.Nth);
            			
            ExtPrb = Solvers.ExteriorSolver( struct(...
                     'Basis'            , Basis.Exterior, ...
                     'Grid'             , Grids.Exterior, ...
                     'CoeffsHandle'     , WaveNumberHandle.Exterior, ... 
                     'CoeffsParams'     , WaveNumberParams.Exterior, ...
                     'ScattererHandle'  , ScattererHandle.Exterior, ...
                     'ScattererParams'  , ScattererParams.Exterior, ...
                     'CollectRhs'       , 0, ... %i.e. no
                     'Extension'        , ExtentionParams.ExtExtension, ...
                     'ExtensionParams'  , [] ...
                     ));

            Grids.Interior = Tools.Grid.CartesianGrid(Boundaries.Interior.x1,Boundaries.Interior.xn,N.Interior.Nx,Boundaries.Interior.y1,Boundaries.Interior.yn,N.Interior.Ny);
            
            IntPrb = Solvers.InteriorHomoSolver( struct(...
                     'Basis'            , Basis.Interior, ...
                     'Grid'             , Grids.Interior, ...
                     'CoeffsHandle'     , WaveNumberHandle.Interior, ... 
                     'CoeffsParams'     , WaveNumberParams.Interior, ...
                     'ScattererHandle'  , ScattererHandle.Interior, ...
                     'ScattererParams'  , ScattererParams.Interior, ...
                     'CollectRhs'       , 1, ... %i.e. yes
                     'Extension'        , ExtentionParams.IntExtension, ...
                     'ExtensionParams'  , [] ...
                     ));

            Grids.Medium = Tools.Grid.CartesianGrid(Boundaries.Medium.x1,Boundaries.Medium.xn,N.Medium.Nx,Boundaries.Medium.y1,Boundaries.Medium.yn,N.Medium.Ny);

            MedPrb = Solvers.InteriorHomoSolver( struct(...
                     'Basis'            , Basis, ...
                     'Grid'             , Grids.Medium, ...
                     'CoeffsHandle'     , WaveNumberHandle.Medium, ... 
                     'CoeffsParams'     , WaveNumberParams.Medium, ...
                     'ScattererHandle'  , ScattererHandle.Medium, ...
                     'ScattererParams'  , ScattererParams.Medium, ...
                     'CollectRhs'       , 1, ... %i.e. yes
                     'Extension'        , Extension, ...
                     'ExtensionParams'  , ExtentionParams ...
                     ));
                          
            if 1
                
                UincParams2  = struct('ScattererType','circle','r',ExtPrb.Scatterer.r);
                
                rhs = zeros(numel(ExtPrb.GridGamma) + numel(MedPrb.GridGamma) + numel(IntPrb.GridGamma),1);
                % rhs(numel(IntPrb.GridGamma)+1:end,1)= -Uinc(UincParams,ExtPrb.Scatterer.th,IncAng,k);
                uinc = Uinc(UincParams2,ExtPrb.Scatterer.th,IncAng.Radians,K{k}.Exterior);
                rhs(numel(IntPrb.GridGamma)+numel(MedPrb.GridGamma)+1:end,1)= ExtPrb.Qcol2(uinc);
                
                [ne,me] = size([ExtPrb.Q{1},ExtPrb.Q{2}]);
                [ni,mi] = size([IntPrb.Q{1},IntPrb.Q{2}]);

                Cn	= [ IntPrb.Q{1}     , IntPrb.Q{2}   , spalloc(ni,me,0)                  ;   
                        MedPrb.Q{1,1}   , MedPrb.Q{1,2} , MedPrb.Q{2,1} , MedPrb.Q{2,2}     ;
                        spalloc(ne,mi,0)                , ExtPrb.Q{1}   , ExtPrb.Q{2} ...
                        ]\rhs ;

                
                %             cn0 = cn(1:2*M+1);
                %             cn1 = cn(2*M+2:end);
                
                Intxi = spalloc(N.Interior.Nx,N.Interior.Ny   ,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = [IntPrb.W{1}(IntPrb.GridGamma,:),IntPrb.W{2}(IntPrb.GridGamma,:)]*Cn(1:Basis.Interior.NBss0 +Basis.Interior.NBss1) ;
                Intu = IntPrb.P_Omega(Intxi);
                
                Medxi = spalloc(N.Medium.Nx,N.Medium.Ny   ,length(IntPrb.GridGamma));
                Medxi(MedPrb.GridGamma) = [MedPrb.W{1,1}(MedPrb.GridGamma,:),MedPrb.W{1,2}(MedPrb.GridGamma,:), ...
                                           MedPrb.W{2,1}(MedPrb.GridGamma,:),MedPrb.W{2,2}(MedPrb.GridGamma,:) ]*Cn;
                Medu = MedPrb.P_Omega(Medxi);
                
                Extxi = spalloc(N.Exterior.Nr,N.Exterior.Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = [ExtPrb.W{1}(ExtPrb.GridGamma,:),ExtPrb.W{2}(ExtPrb.GridGamma,:)]*Cn(Basis.Interior.NBss0 +Basis.Interior.NBss1+1:end) ;%- uinc;
                %Extxi(ExtPrb.GridGamma) = [ExtPrb.W(ExtPrb.GridGamma,:),uinc]*cn;
                
                UincParams3  = struct('ScattererType','circle','r',Grids.Exterior.R);
                uinc = Uinc(UincParams3,Grids.Exterior.Theta,IncAng.Radians,K{k}.Exterior);
                
                
                %tmp = Uinc(UincParams,PlrGrid.Theta,IncAng,k);
%                 ExtuInc = spalloc(PlrGrid.Nx,PlrGrid.Ny,length(ExtPrb.Scatterer.Nm));
%                 ExtuInc(ExtPrb.Scatterer.Outside) = tmp(ExtPrb.Scatterer.Outside);
                Extu = ExtPrb.P_Omega(Extxi,uinc ); 
                
            else
                
              
                cn1 =( ExtPrb.Q1 \ ( -ExtPrb.Q0*Basis.Exterior.cn0 )) ; 
                
                Extxi = spalloc(N.Exterior.Nr,N.Exterior.Nth-1,length(ExtPrb.GridGamma));
                Extxi(ExtPrb.GridGamma) = ExtPrb.W0(ExtPrb.GridGamma,:)*Basis.Exterior.cn0 + ExtPrb.W1(ExtPrb.GridGamma,:)*cn1;
                Extu = ExtPrb.P_Omega(Extxi);
                
                
                cn1 =( IntPrb.Q1 \ ( -IntPrb.Q0*Basis.Interior.cn0 )) ;
                
                Intxi = spalloc(N.Interior.Nx,N.Interior.Ny,length(IntPrb.GridGamma));
                Intxi(IntPrb.GridGamma) = IntPrb.W0(IntPrb.GridGamma,:)*Basis.Interior.cn0 + IntPrb.W1(IntPrb.GridGamma,:)*cn1;
                Intu = IntPrb.P_Omega(Intxi);
                
                
                Cn	= [ MedPrb.Q{1,2}, MedPrb.Q{2,2}] \ ([ -MedPrb.Q{1,1},-MedPrb.Q{2,1}]*[Basis.Interior.cn0; Basis.Exterior.cn0]) ;
                InteriorCn1 = Cn(1:Basis.Interior.NBss1);
                ExteriorCn1 = Cn((1+Basis.Interior.NBss1):end);
                
                Medxi = spalloc(N.Medium.Nx,N.Medium.Ny,length(MedPrb.GridGamma));
                Medxi(MedPrb.GridGamma) = MedPrb.W{1,1}(MedPrb.GridGamma,:)*Basis.Interior.cn0 + MedPrb.W{1,2}(MedPrb.GridGamma,:)*InteriorCn1 ...
                                        + MedPrb.W{2,1}(MedPrb.GridGamma,:)*Basis.Exterior.cn0 + MedPrb.W{2,2}(MedPrb.GridGamma,:)*ExteriorCn1;
                Medu = MedPrb.P_Omega(Medxi);
                
            end
            
            t=toc;
            
            %%%%%%%%%%%%%%
            
           
            
            % % % % % % % % % % % % % % % %
            % Comparison
            % % % % % % % % % % % % % % % %
            if n > 1
                Extu1= spalloc(N.Exterior.Nr,N.Exterior.Nth-1,nnz(Extu0));
                Extu1(ExtPrb.Scatterer.Nm) = Extu0(ExtPrb.Scatterer.Nm);
                
                Medu1= spalloc(N.Medium.Nx,N.Medium.Ny,nnz(Medu0));
                Medu1(MedPrb.Np) = Medu0(MedPrb.Np);
                
                Intu1= spalloc(N.Interior.Nx,N.Interior.Ny,nnz(Intu0));
                Intu1(IntPrb.Np) = Intu0(IntPrb.Np);
                
                
                tmp = Extu(1:2:end,1:2:end)-Extu1(1:2:end,1:2:end);
                ExtErr =norm(tmp(:),inf);

                tmp = Medu(1:2:end,1:2:end)-Medu1(1:2:end,1:2:end);
                MedErr =norm(tmp(:),inf);
                
                tmp = Intu(1:2:end,1:2:end)-Intu1(1:2:end,1:2:end);
                IntErr =norm(tmp(:),inf);
                
                TotErr = max([IntErr,MedErr,ExtErr]);
                
                fprintf('N.Int=%-5dx%-5d,\t N.Med=%-5dx%-5d\t N.Ext=%-5dx%-5d,\t  TotErr: %d|%-5.2f\t IntErr: %d|%-5.2f\t MedErr: %d|%-5.2f\t ExtErr: %d|%-5.2f\t time=%d\n', ...
                                    N.Interior.Nx,N.Interior.Ny,N.Medium.Nx,N.Medium.Ny,N.Exterior.Nr,N.Exterior.Nth,TotErr,log2(TotErrPre/TotErr),IntErr,log2(IntErrPre/IntErr),MedErr,log2(MedErrPre/MedErr),ExtErr,log2(ExtErrPre/ExtErr),t);
                                
                TotErrPre = TotErr; IntErrPre=IntErr; MedErrPre=MedErr; ExtErrPre=ExtErr;
            end
                        
            Extu0=spalloc(N.Exterior.Nr*2-1,N.Exterior.Nth*2-2,nnz(Extu));
            Extu0(1:2:end,1:2:end)=Extu;

            Medu0=spalloc(N.Medium.Nx*2-1,N.Medium.Ny*2-1,nnz(Intu));
            Medu0(1:2:end,1:2:end)=Medu;
            
            Intu0=spalloc(N.Interior.Nx*2-1,N.Interior.Ny*2-1,nnz(Intu));
            Intu0(1:2:end,1:2:end)=Intu;

            if want2plot || n==nmax
                
                
                [XExt,YExt,tExtu] = Tools.Common.MaskU(Grids.Exterior,ExtPrb.Scatterer.Mp,Extu);
                [XInt,YInt,tIntu] = Tools.Common.MaskU(Grids.Interior,IntPrb.Scatterer.Mm,Intu);
                [XMed,YMed,tMedu] = Tools.Common.MaskU(Grids.Medium  ,MedPrb.Scatterer.Mm,Medu);
                
                figure; Tools.Common.mPcolor([XInt,XMed,XExt],[YInt,YMed,YExt],full(abs([tIntu,tMedu,tExtu])),'total field, abs')
                %figure; Tools.Common.mPcolor([XInt,XMed,XExt],[YInt,YMed,YExt],full(real([tIntu,tMedu,tExtu])),'total field, real')
                %figure; Tools.Common.mPcolor([XInt,XMed,XExt],[YInt,YMed,YExt],full(imag([tIntu,tMedu,tExtu])),'total field, imag')
                
                if 1 % debug, draw wave number
                    [XExt,YExt,KExt] = Tools.Common.MaskU(Grids.Exterior,ExtPrb.Scatterer.Mp,ExtPrb.k*ones(size(Grids.Exterior.R)));
                    KInt=WaveNumberHandle.Interior(Grids.Interior,WaveNumberParams.Interior);
                    [XInt,YInt,KInt] = Tools.Common.MaskU(Grids.Interior,IntPrb.Scatterer.Mm,KInt.k);
                    KMed=WaveNumberHandle.Medium(Grids.Medium,WaveNumberParams.Medium);
                    [XMed,YMed,KMed] = Tools.Common.MaskU(Grids.Medium,MedPrb.Scatterer.Mm,KMed.k);
                    figure
                    mesh([XInt,XMed,XExt],[YInt,YMed,YExt],[KInt,KMed,KExt])
                    
                    T = sprintf('k: Int:,%d, Med: %d, Ext %d \n',  K{k}.Interior , K{k}.Medium, K{k}.Exterior);
                    title(T);
                end
                
                %               

                %tmp=abs([tIntu,tMedu,tExtu]);
                %fprintf('min=%d \t max=%d\n', full(min(tmp(:))), full(max(tmp(:))));
                
                % saveas(figure(1),'totfldkin15kout30abs.jpg','jpg')
                % saveas(figure(2),'totfldkin15kout30real.jpg','jpg')
                % saveas(figure(3),'totfldkin15kout30imag.jpg','jpg')
            end
            
            
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



