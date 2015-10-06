function RunNestedHomo
    % Semyon Method, Homo HLM, Cart coord, 
    


x1=-2.2;xn=2.2;
y1=-2.2;yn=2.2;

R0=1;
R1=2;
NHR=1.6;


ScatType = 'circle'; %'ellipse';% 'StarShapedScatterer'; %'ellipse' or 'circle' or 'StarShapedScatterer'
BType = 'Fourier'; % 'Fourier' or 'Chebyshev'
ChebyshevRange = struct('a',-pi,'b',pi);%don't change it

if strcmpi(ScatType,'circle')
    ExParams{1}  = struct('ScattererType','circle', 'r', R0);
    ExParams{2}  = struct('ScattererType','circle', 'r', R1);

end

%fprintf('Method:%s,\t Ellipse: a=%d; \t b=%d \n',ScatType,a,b);
fprintf('Method:RunInteriorHomo-%s,\t  \n',ScatType);
fprintf('Grid: x1=%f, xn=%f, y1=%f, yn=%f \n , R0=%f, R1=%f \n',x1,xn,y1,yn,R0,R1);

for k =1%[1,5,10,15,20,25]
    
    ErrPre = 0; ErrXi2Pre =0; ErrXi1Pre =0;
    
    f   =@(th,n) Exact  (th,k,ExParams{n});
    dfdn=@(th,n) drExact(th,k,ExParams{n});
	
	if strcmpi(BType,'Chebyshev')
		Basis = struct( 'Interior', Tools.Basis.ChebyshevBasis.BasisHelper(@(th) f(th,1),@(th) dfdn(th,1),ChebyshevRange), ...
                        'Exterior', Tools.Basis.ChebyshevBasis.BasisHelper(@(th) f(th,2),@(th) dfdn(th,2),ChebyshevRange) ...
                       );
    elseif strcmpi(BType,'Fourier')
		Basis = struct( 'Interior', Tools.Basis.FourierBasis.BasisHelper(@(th) f(th,1),@(th) dfdn(th,1)), ...
                        'Exterior', Tools.Basis.FourierBasis.BasisHelper(@(th) f(th,2),@(th) dfdn(th,2)) ...
                       );
	end
fprintf('NBss: %d,%d,%d,%d\n',Basis.Interior.NBss0,Basis.Interior.NBss1,Basis.Exterior.NBss0,Basis.Exterior.NBss1);

for n=1:4 %run different grids
tic
    p=4;%3;
	Nx=2.^(n+p)+1;	Ny=2.^(n+p)+1;
    
    if strcmpi(ScatType,'circle')
        ScattererHandle  = @Tools.Scatterer.NestedPolarScatterer;
        ScattererParams  = struct('r0',R0,'r1',R1, 'Stencil', 9);
        Extension        = @Tools.Extensions.NestedExtension;
        ExtentionParams  = struct('IntExtension',@Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension,'ExtExtension',@Tools.Extensions.EBPolarHomoHelmholtz5OrderExtension);
    end
    
    IntPrb =  Solvers.InteriorHomoSolver( struct(...
        'Basis'          , Basis, ...
        'Grid'           , Tools.Grid.CartesianGrid(x1,xn,Nx,y1,yn,Ny), ...
        'CoeffsHandle'   , @Tools.Coeffs.ConstantWaveNumber, ... 
        'CoeffsParams'   , struct('k',k,'r0',NHR), ...
        'ScattererHandle', ScattererHandle, ...
        'ScattererParams', ScattererParams, ...
        'CollectRhs'     , 1, ... %i.e. yes
        'Extension'      , Extension, ...
        'ExtensionParams', ExtentionParams ...
        ));
    
    InteriorCn1	= ( IntPrb.Q{1,2} \ ( -IntPrb.Q{1,1}*Basis.Interior.cn0)) ;
    ExteriorCn1 = ( IntPrb.Q{2,2} \ ( -IntPrb.Q{2,1}*Basis.Exterior.cn0)) ;
    

    %Cn	= [ IntPrb.Q{1,2}, IntPrb.Q{2,2}] \ ([ -IntPrb.Q{1,1},-IntPrb.Q{2,1}]*[Basis.Interior.cn0, Basis.Exterior.cn0]) ;
    %InteriorCn1 = Cn(1:Basis.Interior.NBss0);
    %ExteriorCn1 = Cn(1+Basis.Interior.NBss0:end);
    
    xi = spalloc(Nx,Ny,length(IntPrb.GridGamma));
    xi(IntPrb.GridGamma{1}) = IntPrb.W{1,1}.W(IntPrb.GridGamma{1},:)*Basis.Interior.cn0 + IntPrb.W{1,2}.W(IntPrb.GridGamma{1},:)*InteriorCn1;
    xi(IntPrb.GridGamma{2}) = IntPrb.W{2,1}.W(IntPrb.GridGamma{2},:)*Basis.Exterior.cn0 + IntPrb.W{2,2}.W(IntPrb.GridGamma{2},:)*ExteriorCn1;

    
    if strcmpi(ScatType,'circle')
        ExParams2 =ExParams;
        ExParams2{1}.r = IntPrb.Scatterer.r{1};
        ExParams2{2}.r = IntPrb.Scatterer.r{2};
        xiex1 = Exact(IntPrb.Scatterer.th{1},k,ExParams2{1});
        xiex2 = Exact(IntPrb.Scatterer.th{2},k,ExParams2{2});
    end
    
    ErrXi1 =norm(xiex1 -xi(IntPrb.GridGamma{1}),inf);
    ErrXi2 =norm(xiex2 -xi(IntPrb.GridGamma{2}),inf);
    %xi(IntPrb.GridGamma{1})=xiex1;
    %xi(IntPrb.GridGamma{2})=xiex2;
    u = IntPrb.P_Omega(xi);
    
    t1=toc;
    
    % % % % % % % % % % % % % % % %
    % Comparison
    % % % % % % % % % % % % % % % %
    
    
    exact = zeros(Nx,Ny);   
    %exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.R(IntPrb.Scatterer.Np),IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k);
    if strcmpi(ScatType,'ellipse')
        ExParams3  = struct('ScattererType','ellipse','eta',IntPrb.Scatterer.Eta(IntPrb.Scatterer.Np),'FocalDistance',FocalDist,'ExpansionType',25, 'Stencil', 9);
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Phi(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
        %TBD;
    elseif strcmpi(ScatType,'circle')
       % ExParams3 =ExParams;
        ExParams3.r = IntPrb.Scatterer.R(IntPrb.Scatterer.Np);
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);%(IntPrb.Scatterer.r,IntPrb.Scatterer.th,k);
     elseif strcmpi(ScatType,'StarShapedScatterer')
        
        ExParams3 = struct('ScattererType','StarShapedScatterer','r',IntPrb.Scatterer.R(IntPrb.Scatterer.Np), 'Stencil', 9);
        
        exact(IntPrb.Scatterer.Np) = Exact(IntPrb.Scatterer.Th(IntPrb.Scatterer.Np),k,ExParams3);        
    
    end
   
    
    
    %t2=toc;
    
    ErrTot =norm(exact(:)-u(:),inf);
    fprintf('k=%d,N=%-4dx%-4d\t ErrXi=%d|%d\t rate=%-5.2f|%-5.2f ErrTot=%d\t rate=%-5.2f timeA=%d\n',...
        k, Nx,Ny,ErrXi1,ErrXi2,log2(ErrXi1Pre/ErrXi1),log2(ErrXi2Pre/ErrXi2),ErrTot,log2(ErrPre/ErrTot),t1);
    ErrPre = ErrTot;
    ErrXi1Pre = ErrXi1;
    ErrXi2Pre = ErrXi2;
    
   % figure, plot(x,abs(u(fix(Ny/2),:)),'r+-',x,abs(exact(fix(Ny/2),:)),'bo');
%     saveas(gcf,['inhomo_k=' num2str(k) 'N=' num2str(Nx) '.jpg'],'jpg');
    % plot(x,abs(u(:,fix(Ny/2))),'r+-',x,abs(exact(:,fix(Ny/2))),'bo');

end

fprintf('\n');

%Linf=log2(etinf(1:end-1)./etinf(2:end))
%Lbinf=log2(ebinf(1:end-1)./ebinf(2:end))
end
end



function g = Exact(th,k,Params)%(R,th,k)
% if strcmpi(Params.ScattererType,'ellipse')
%     x = Params.FocalDistance * cosh(Params.eta) .* cos(th);
%     %          y = Params.FocalDistance * sinh(Params.eta) .* sin(th);
% elseif strcmpi(Params.ScattererType,'circle')
    x = Params.r .* cos(th);
    %         y = Params.r .* sin(th);
% elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
%     try
%         x = Params.Parameterization.XHandle.Derivatives(th);
%     catch
%         x = Params.r.*cos(th);
%     end
% end

%g = exp(1i.*k.*R.*cos(th));
g = exp(1i.*k.*x);
end


function dne = drExact(th,k,Params)%(R,th,k)
if strcmpi(Params.ScattererType,'ellipse')
    dx = Params.FocalDistance * sinh(Params.eta) .* cos(th);
    %         dy = Params.FocalDistance * cosh(Params.eta) .* sin(th);
elseif strcmpi(Params.ScattererType,'circle')
    dx = cos(th);
    %         dy = sin(th);
elseif strcmpi(Params.ScattererType,'StarShapedScatterer')
    [x,dx] = Params.Parameterization.XHandle.Derivatives(th);
    [y,dy] = Params.Parameterization.YHandle.Derivatives(th);
end

e = Exact(th,k,Params);%(R,th,k);
%dne = 1i.*k.*cos(th).*e;
dne = 1i.*k.*dx.*e;


if strcmpi(Params.ScattererType,'StarShapedScatterer')
   h = sqrt(dx.^2 + dy.^2);
   % dne = dne./h;
	
	%fx dy + fy (-dx)
	
	dne  = (1i.*k.*e).*dy./h; %=fx dy, here fy=0
	
end


end
