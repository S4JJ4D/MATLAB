function [] = quadPlot_1D(simDataset)

%% 3D Model Demo
% This is short demo that loads and renders a 3D model of a human femur. It
% showcases some of MATLAB's advanced graphics features, including lighting and
% specular reflectance.

% Copyright 2011 The MathWorks, Inc.


%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
% fv = stlread('Assembly.stl');
load('STL_Data');

%% Initialization
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.

%%Creating figure with customized specifications
close all;
fg = figure(1);
fg.Units = 'Normalized';
fg.Position = [0.1740 0.2142 0.5786 0.6050];
ax = axes(fg);
quadPatch = patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
    'EdgeColor',       'none',        ...
    'FaceLighting',    'gouraud',     ...
    'AmbientStrength', 0.15);
hold on;
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
% view([135 35]);
view(0,0);
grid on
axis([-1 1 -1 1 0 2]);
%%Labeling coordinate frame Axes
xlabel('X','FontWeight','bold','FontSize',12);
ylabel('Y','FontWeight','bold','FontSize',12);
zlabel('Z','FontWeight','bold','FontSize',12);

%%Plotting dashed lines aligned with global axes
plot3([-100;100],[0;0],[0;0],'r-.');
plot3([0;0],[100;-100],[0;0],'r-.');
plot3([0;0],[0;0],[100;-100],'r-.');
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[.3;0;0],[0;.3;0],[0;0;.3],'LineWidth',1.5,'Color','blue');
text(ax,'String','z_0','Position',[0,0,.3],'FontWeight','bold');
text(ax,'String','y_0','Position',[0,.3,0],'FontWeight','bold');
text(ax,'String','x_0','Position',[.3,0,0],'FontWeight','bold');

%%Plotting dashed lines aligned with body-fixed axes
arrowLength = .5;
h = quiver3(zeros(3,1),zeros(3,1),zeros(3,1),...
    [arrowLength;0;0],[0;arrowLength;0],[0;0;arrowLength],...
    'LineWidth',1.5,'Color','black');
xyPlot = plot3([0 0 0],[0 0 0],[0 0 0],':ok');

%%Creaing a transofrm object and assigning its children:
trObject = hgtransform;
h.Parent = trObject;
quadPatch.Parent = trObject;

Thrust_group = hggroup;
Thrust_group.Parent = trObject;


%%Labeling body-fixed axes:
text(ax,'String','z_b','Position',[0,0,arrowLength],'Parent',trObject','FontWeight','bold','FontSize',12);
text(ax,'String','y_b','Position',[0,arrowLength,0],'Parent',trObject','FontWeight','bold','FontSize',12);
text(ax,'String','x_b','Position',[arrowLength,0,0],'Parent',trObject','FontWeight','bold','FontSize',12);

%%Labeling quadrotor arms:
wingSpan = .175;
textWingspan = .2;
offset = .05;
text(ax,'String','A','Position',[0,-textWingspan,offset],'Parent',trObject','FontWeight','bold','FontSize',12);
text(ax,'String','B','Position',[0, textWingspan,offset],'Parent',trObject','FontWeight','bold','FontSize',12);
text(ax,'String','C','Position',[textWingspan,0,offset],'Parent',trObject','FontWeight','bold','FontSize',12);
text(ax,'String','D','Position',[-textWingspan,0,offset],'Parent',trObject','FontWeight','bold','FontSize',12);

%Visualizing relative forces produced by motors:
F_A = quiver3(0,-wingSpan,0,0,0,.2,...
    'LineWidth',2.5,'Color','red','Parent',Thrust_group);
F_B = quiver3(0, wingSpan,0,0,0,.2,...
    'LineWidth',2.5,'Color','red','Parent',Thrust_group);
F_C = quiver3(wingSpan,0,0,0,0,.2,...
    'LineWidth',2.5,'Color','red','Parent',Thrust_group);
F_D = quiver3(-wingSpan,0,0,0,0,.2,...
    'LineWidth',2.5,'Color','red','Parent',Thrust_group);


%%Creating text annotation:
Textannot = annotation(fg,'textbox',...
    [0.45 0.93 0.11 0.044],...
    'String','Z = 0,    time = 0',...
    'LineStyle','none',...
    'FontSize',15,...
    'FitBoxToText','on','FontWeight','bold');


time = simDataset{1}.Values.Time;
freq = 1/median(diff(time));

z =  simDataset{1}.Values.Data;
thrust_cmd = simDataset{2}.Values.Data; thrust_cmd = thrust_cmd/(2.5*max(thrust_cmd));

figure(fg);
pause(1);

%% Iterative Plot

t=0; %time initialization
u = @(p) p>0;
n = 1;
tic;    %start time measuring
% I added the "ishandle" so the program will end in case u closed the figure
while (n <= numel(time))
    M = makehgtform('translate',[0,0,z(n)]);
    trObject.Matrix = M;
%     xyPlot.XData = [0,x(n),x(n)];
%     xyPlot.YData = [0,y(n),y(n)];
    xyPlot.ZData = [0,0,z(n)];
    drawnow;  %updates the display
    axis([-1,1,-1,1,u(-z(n))*z(n),u(z(n))*z(n)+1]);
    set(Thrust_group.Children(:),'WData',thrust_cmd(n));
    Textannot.String = sprintf('Z = %.2f,   time = %.2f : %d',z(n),t,time(end));
    t = toc; %measuring current time
    n = round(t*freq)+1;    
end
Textannot.String = sprintf('Z = %.2f,   time = %.2f : %d',z(numel(time)),time(end),time(end));
trObject.Matrix = M;

end