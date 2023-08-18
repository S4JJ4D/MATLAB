% RUN ME!

clc;
clear;
close all;

if exist('fg', 'var')
    if ishandle(fg)
        delete(fg);
    end
end

fg = figure(...
    'Units', 'normalized', ...
    'Position', [0.6453 0.3417 0.3464 0.5444]);
ax = axes;
hold(ax, 'on');
grid on;
box on;
% xline(0, 'k--');
% yline(0, 'k--');
axis([-1 11 -8 2]);
ax.DataAspectRatio = [1 1 1];


%%
offset = 4;
f = [1, 2, 3, 4];
x = [0, 1.5, 1.5, 0] + offset;
y = [-1, -1, 1, 1];

types_list = ["simplified", "detailed", "axled"];
numTypes = numel(types_list);

boxes = gobjects(1, numTypes);

for i=1:numTypes
    boxes(i) = patch(ax, 'XData', x, 'YData', y - 3*(i-1), ...
    'FaceColor', 'k', 'FaceAlpha', .5, 'LineWidth', 1.5);

    springs(i) = DynamicSpring("spring"+i, 'Radius', .25, 'Pitch', .03, 'Turns', 18,...
    'Axes', ax, 'VisualForm', types_list(i));
    springs(i).plotting_options.FrontEyeInnerColor = 'r';

    texts(i) = text(8, -3*(i-1), types_list(i), ...
        'FontName', 'Source Code Pro', 'FontWeight', 'bold', ...
        'FontSize', 11); 
end


%% Simulation and Animation

m = .05;
b = .05;
k = 2.5;
sys = tf(1, [m b k]);
time = 0:.01:5;
x_out = impulse(sys, time);

Qx = x_out.' + offset;
Qy = zeros(1, numel(time));

n=1;
freq = 1/median(diff(time));


for i=1:numTypes
    springs(i).PlotSpring(vecnorm([Qx(n), Qy(n)]), 'Configuration', [0, -3*(i-1), atan2(Qy(n), Qx(n))]);
end

pause(.5);

t = 0; %time initialization
n = 1;
tic;    %start time measuring
while (n <= numel(time))

    for i=1:numTypes
        set(boxes(i), 'XData', [0, 1.25, 1.25, 0] + Qx(n));
        springs(i).PlotSpring(vecnorm([Qx(n), Qy(n)]), 'Configuration', [0, -3*(i-1), atan2(Qy(n), Qx(n))]);
    end

    drawnow;
    t = toc; %measuring current time
    n = round(t*freq)+1;
end
