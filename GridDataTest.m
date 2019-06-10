%test for Extrapolation to gamma, used in Descrete Source...

errp = 0;
R0=1;
f = @(x,y) sin(x).*cos(y);

for n=3:9
    N = 2^n;
    Grid = Tools.Grid.CartesianGrid(-2,2,N,-2,2,N);
    
    if 1
        PS=Tools.Scatterer.PolarScatterer(Grid,struct('r0',R0, 'Stencil', 5));
        
        Params.Src = f(Grid.X,Grid.Y);
        Params.Np = PS.Np;
        Params.Mp = PS.Mp;
        Params.GridGamma = PS.GridGamma;
        Params.Grid = Grid;
        S = Tools.Source.NavierStokesDescreteSrc(PS.TheScatterer(),[],[],Params);
        Nu = S.Derivatives();
        
        th = PS.th;
    else
        th=linspace(0,2*pi,N);
                
        U = f(Grid.X,Grid.Y);
        
        msk = Grid.R<R0;
        
        Uu=U(msk);
        X = Grid.X(msk);
        Y = Grid.Y(msk);
        
        %  Nu = interp2(X,Y,Uu,nx,ny,'makima');%'cubic');
        
        Nu = griddata(X,Y,Uu,nx,ny,'v4');%'cubic');
        %        F = scatteredInterpolant(X,Y,Uu,'linear');
        %        Nu=F(nx,ny);
        
        %     R=Grid.R(Grid.R<1);
        %     Theta=Grid.Theta(Grid.R<1);
        %     F = scatteredInterpolant(R,Theta,Uu,'natural','linear');
        %     Nu=F(ones(size(th))*R0,th);
        
        %     spline = csapi({X,Y},Uu);
        
    end
    
    nx = R0*cos(th);
    ny = R0*sin(th);
    
    err = norm( Nu - f(nx,ny),inf);
    
    fprintf('%f\n',log2(errp/err))
    
    errp = err;
    
end