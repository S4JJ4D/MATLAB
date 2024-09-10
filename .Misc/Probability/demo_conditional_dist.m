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
ax = axes;

% 1. Define the distribution object:
N1 = makedist('Normal','mu',0,'sigma',2.5);
N2 = makedist('Normal','mu',9,'sigma',4);
N3 = makedist('Normal','mu',15,'sigma',1.5);

% plot PDF:
x = -10:.01:25;
y1 = pdf(N1, x);
y2 = pdf(N2, x);
y3 = pdf(N3, x);
% plot(x, y1, 'LineWidth', 1);
hold on;
% plot(x, y2, 'LineWidth', 1);
% plot(x, y3, 'LineWidth', 1);


f1 = @(x) 1/(sqrt(2*pi)*2.5) * exp(-1/2 * ((x-0)/2.5).^2);
f2 = @(x) 1/(sqrt(2*pi)*4) * exp(-1/2 * ((x-9)/4).^2);

ft = @(x) .5*f1(x) + .5*f2(x);
area_total = integral(ft,-10,25);

b = 4;
area_b = integral(ft, b, 25);

ft2 = @(x) ft(x)/area_b;
xx = b:.01:25;
yt2 = ft2(xx);
area(xx, yt2, 'FaceColor', '#DAD9EB', 'EdgeColor', '#4A457F', 'FaceAlpha', .5, ...
    'LineWidth', 1);
line([b, b], [0, ft2(b)], 'Color', '#4A457F', 'LineWidth', 1);

yt = .5*y1 + .5*y2;
area(x, yt, 'FaceColor', '#F6E1BE', 'EdgeColor', '#E5A73C', 'FaceAlpha', .5, ...
    'LineWidth', 1);

X_population = x;
pdf = yt;

set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Box', 'on', ...
    'FontName', 'Source Code Pro');

t1 = text(11.5, .09, '$$\leftarrow f_{X|A}(x|A)$$', 'Interpreter', 'latex', 'FontSize', 11);
t2 = text(-6.1, 0.07, '$$f_{X}(x)\rightarrow$$', 'Interpreter', 'latex', 'FontSize', 11);


