function p = plotcube(X,Y,Z)
% PLOTCUBE - Display a 3D-cube in the current axes

x = X/2;
y = Y/2;
z = Z/2;

vert = [x,y,-z;-x,y,-z;-x,-y,-z;x,-y,-z;x,y,z;-x,y,z;-x,-y,z;x,-y,z];
fac = [1 2 3 4;1 2 6 5;1 4 8 5;2 3 7 6;3 4 8 7;5 6 7 8];
p = patch('Vertices',vert,'Faces',fac,'FaceColor','red');
view(3);

% 
% h = hggroup;
% for i=1:6, L(i).Parent = h; end

view(3);
