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


clear;
close all;

fg = figure;


x = 0:.01:1;
y = x.^2;

area(x, y, 'FaceColor', '#DAD9EB', 'EdgeColor', '#4A457F');
hold on;
ax = gca;

set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XAxisLocation', 'origin', 'YAxisLocation', 'origin', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');

xlabel('x');
ylabel('y');

t = text(.43, .4, '$$y=x^2 \rightarrow$$', 'Interpreter', 'latex');
vline = plot([1, 1], [0, 1], 'Color', '#4A457F');
axis equal;
axis([-0.4487    1.0855   -0.1507    1.0593]);


% best mse predictor

x = 0:.01:2;
y_hat_1 = (x.^2)/2;
y_hat_2 = 2/3 * x - 1/5;


p1 = plot(x, y_hat_1, 'Color', [0.7412, 0.0588, 0.0588],...
    'LineWidth', 1, 'DisplayName', 'Best MSE Predictor');
p2 = plot(x, y_hat_2, 'k-', 'LineWidth', 1, 'DisplayName', 'Best Linear Predictor');

legend([p1,p2]);




