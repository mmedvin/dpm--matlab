Exact=Tools.Exact.NSExact1;

Nt=1000;

Stencil = 5;
ScattererParams = struct('r0',0.9, 'Stencil', Stencil);

ExParams.r0 = 1; % = struct(','r',Grid.R);
Eparam      = @(r) struct('r0',ExParams.r0,'r',r);

%Ft = @smoothstep;
Ft = @cossin;

if 0 
    fprintf('Starting Test 1\n')    

    preerr=0;
    
    for  n=3:11
        Grid = Tools.Grid.CartesianGrid(-1.2,1.2,2^n,-1.2,1.2,2^n);
        ExParams.r = Grid.R;
        
        %O=Exact.Omega(Grid.Theta, ExParams);
        %LO=Tools.Exact.NSExact1.LaplacianOmega(Grid.Theta, ExParams);
        %[ft,df]=smoothstep(1.3);
        
        RN=10;
        %f=df*O-LO*ft/RN;
        %F=Exact.NSSource(Grid.Theta, ExParams,RN,@smoothstep,1.3,0);
        %mesh(Grid.X,Grid.Y,real(f-F))
        
        S = Tools.Scatterer.PolarScatterer(Grid,ScattererParams);
        
        T=zeros(Grid.Nx,Grid.Ny);
        Msk = S.Np;
        T(Msk)=Exact.Test(Grid.Theta(Msk), Eparam(Grid.R(Msk)),RN,Grid.dx,@smoothstep,10,0);
        
        
        %T(Grid.R>1.1)=0;
        %mesh(Grid.X,Grid.Y,abs(T))
        
        err = max(T(:));
        fprintf('N=%d, e=%d r=%d\n',2^n,err,log2(preerr/err))
        preerr=err;
    end
    
fprintf('End of Test 1\n')    
end


if 0
    fprintf('Testing Laplacian (numerical vs exact) \n')    

    preerr=0;
    
    for  n=3:11
        Grid = Tools.Grid.CartesianGrid(-1.2,1.2,2^n,-1.2,1.2,2^n);
        ExParams.r = Grid.R;
        
        
        Der = Tools.Common.SecondDerivative(Grid.Nx,Grid.Ny,Grid.dx,Grid.dy);
        O=Exact.Omega(Grid.Theta, ExParams);
        LOe=Exact.LaplacianOmega(Grid.Theta, ExParams);
        LOa=Der.CartesianLaplacian(O);
        
        S = Tools.Scatterer.PolarScatterer(Grid,ScattererParams);
        
        Msk = S.Np;
        
        diff = zeros(Grid.Nx,Grid.Ny);
        diff(Msk) = LOe(Msk) - LOa(Msk);
        
        %mesh(Grid.X,Grid.Y,abs(diff))
        
        err=norm(diff(:),inf);
        
        fprintf('N=%d, e=%d r=%d\n',2^n,err,log2(preerr/err))
        preerr=err;
    end
       fprintf('End of Laplacian test\n')    
 
end

if 1
        fprintf('Starting Test 2\n')    

    preerr=0;
    
    for  n=3:11
        Grid = Tools.Grid.CartesianGrid(-1.2,1.2,2^n,-1.2,1.2,2^n);
        ExParams.r = Grid.R;
                
        RN=100;
                
        S = Tools.Scatterer.PolarScatterer(Grid,ScattererParams);
        
        T=zeros(Grid.Nx,Grid.Ny);
        Msk = 1:numel(T);%S.Np;
        T(Msk)=Exact.Test2(Grid.Theta(Msk), Eparam(Grid.R(Msk)),RN,Grid.dx,Ft,1000,0);
                
        err = norm(T(:),inf);
        fprintf('N=%d, e=%d r=%d\n',2^n,err,log2(preerr/err))
        preerr=err;
    end
        fprintf('End of Test 2\n')    

end
