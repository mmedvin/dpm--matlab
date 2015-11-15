function [ X,Y,U ] = MaskU( Grid,Mask, U)
if isa(Grid, 'Tools.Grid.PolarGrids')
    R=Grid.R;
    Th=Grid.Theta;
    R(:,Grid.Ny)=R(:,Grid.Ny-1);
    Th(:,Grid.Ny)=2*pi;

    X = R .* cos(Th);
    Y = R .* sin(Th);
    if numel(U)==1
        U=U*ones(size(X));
    else
        U(:,Grid.Ny)=U(:,1);
    end
elseif isa (Grid, 'Tools.Grid.CartesianGrid')
    X=Grid.X;
    Y=Grid.Y;
end

X(Mask)=NaN;
Y(Mask)=NaN;

U(Mask)=NaN;
end

