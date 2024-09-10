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

xx = [2, 6, 6, 0, 0, 2];
yy = [0, 0, 6, 6, 2, 0];

close all;
figure;
ax = axes;

pp = patch('XData', xx, 'YData', yy, ...
    'FaceColor', '#DAD9EB', 'EdgeColor', '#4A457F', 'FaceAlpha', .5, 'LineWidth', 1);

set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XAxisLocation', 'origin', 'YAxisLocation', 'origin', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');


tx = xlabel('X', 'Interpreter', 'latex');
ty = ylabel('Y', 'Interpreter', 'latex');

% plot line:
hold on;

ll = plot([-2 4],[4, -2], 'LineStyle', '--', 'Color', '#4A457F');

axis([-4.1754    5.9678   -2.0187    5.9813]);
axis equal;
% intersection region

xx = [1, 6, 6, 1];
yy = [1, 1, 6, 6];

A_region = patch('XData', xx, 'YData', yy, ...
    'FaceColor', '#F6E1BE', 'EdgeColor', '#E5A73C', 'FaceAlpha', .3, 'LineWidth', 1);

A_label = text(3.5,3.5, 'A', 'Interpreter', 'latex');

A_hatch = hatchfill2(A_region,'single','HatchAngle',45);
% hatchfill2(A_region,'cross','HatchAngle',45,'HatchDensity',20,'HatchColor','b','HatchLineWidth',2);


axis equal;

%%
pause(2);

% 
delete([A_region, A_label, A_hatch]);
xx = [1, 1, 6, 6, 0, 0];
yy = [1, 0, 0, 6, 6, 1];

B_region = patch('XData', xx, 'YData', yy, ...
    'FaceColor', '#F6E1BE', 'EdgeColor', '#E5A73C', 'FaceAlpha', .05, 'LineWidth', 1);

text(3,3, 'B', 'Interpreter', 'latex');
B_hatched = hatchfill2(B_region,'single','HatchAngle',45, 'HatchLineWidth',.5);

axis equal;


% #F6E1BE
% #F9EBD3


% #E5A73C
% #E3A537
