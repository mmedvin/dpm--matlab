function [ X,Y,U ] = MaskU( Grid,Mask, U)
if isa(Grid, 'Tools.Grid.PolarGrids')
    X = Grid.R .* cos(Grid.Theta);
    Y = Grid.R .* sin(Grid.Theta);
    if numel(U)==1
        U=U*ones(size(X));
    end
elseif isa (Grid, 'Tools.Grid.CartesianGrid')
    X=Grid.X;
    Y=Grid.Y;
end

X(Mask)=NaN;
Y(Mask)=NaN;

U(Mask)=NaN;
end

