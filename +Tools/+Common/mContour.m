function mContour(X,Y,Z,Title)

contourf(X,Y,Z)

title(Title);
axis equal
axis off;
view(2);
shading flat;
grid off;
H=gca;
set(H,'Color','none');
set(H,'Visible','off');
colorbar
end