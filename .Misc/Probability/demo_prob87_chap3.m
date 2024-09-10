% from Olofsson.

% Mathematica Style:
% 'FontName', 'Source Code Pro'
% Area Colors: 'FaceColor', '#DAD9EB', 'FaceColor', '#F6E1BE'
% Edge Colors: 'EdgeColor', '#4A457F', 'EdgeColor', '#E5A73C'
% Axes Style:
% set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
%     'XAxisLocation', 'origin', 'YAxisLocation', 'origin', ...
%     'Box', 'on', ...
%     'FontName', 'Source Code Pro');
%
%% 
clear;
close all;

xx = [0 0 0 0 1 1 1 2 2 3];
yy = [0 1 2 3 0 1 2 0 1 0];

x_ = [1 1 2];
y_ = [1 2 1];

close all;
figure1 = figure;
ax = axes;

scatter(xx, yy, 'filled');
hold on;
plot(x_, y_, 'o', 'MarkerSize', 13, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
xline(0, 'k-', 'Alpha', .2);
yline(0, 'k-', 'Alpha', .2);

set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');

axis equal;
grid on;

% Uncomment the following line to preserve the X-limits of the axes
xlim(ax,[-0.853167057060849 4.71580226809253]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(ax,[-0.492737052078199 3.8995629479218]);

xlabel('H');
ylabel('S');

% Create arrow
annotation(figure1,'arrow',[0.608571428571428 0.551428571428571],...
    [0.386666666666667 0.386666666666667],'HeadStyle','deltoid');

% Create arrow
annotation(figure1,'arrow',[0.3875 0.3875],[0.656142857142857 0.605714285714286],...
    'HeadStyle','deltoid');

% Create arrow
annotation(figure1,'arrow',[0.4275 0.4025],[0.445238095238096 0.406190476190477],...
    'HeadStyle','deltoid');

% Create textbox
annotation(figure1,'textbox',...
    [0.434214285714285 0.457142857142859 0.172928571428572 0.148095238095243],...
    'VerticalAlignment','middle',...
    'String','$$\frac{\pmatrix{13 \cr 1} \pmatrix{13 \cr 1} \pmatrix{26 \cr 1}}{\pmatrix{52 \cr 3}} $$',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',8,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.326357142857142 0.668571428571431 0.1265 0.149047619047624],...
    'VerticalAlignment','middle',...
    'String','$$\frac{\pmatrix{13 \cr 1} \pmatrix{13 \cr 2}}{\pmatrix{52 \cr 3}} $$',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',8,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.618499999999998 0.311428571428574 0.1265 0.149047619047625],...
    'VerticalAlignment','middle',...
    'String','$$\frac{\pmatrix{13 \cr 1} \pmatrix{13 \cr 2}}{\pmatrix{52 \cr 3}} $$',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',8,...
    'FitBoxToText','off');


