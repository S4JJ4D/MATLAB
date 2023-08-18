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
axis([-1 10 -5 5]);
ax.DataAspectRatio = [1 1 1];


%%
offset = 4;
f = [1, 2, 3, 4];
x = [0, 1.5, 1.5, 0] + offset;
y = [-1, -1, 1, 1];

box_patch = patch(ax, 'XData', x, 'YData', y, ...
    'FaceColor', 'k', 'FaceAlpha', .5, 'LineWidth', 1.5);
hold on;


sp1 = DynamicSpring('spring1', 'Radius', .25, 'Pitch', .03, 'Turns', 18,...
    'Axes', ax, 'VisualForm', 'with_damper');

sp1.plotting_options.FrontEyeInnerColor = 'r';

%% Run Simulation
% loads pre-computed trajectory variables: Qx, Qy, time
% load('motionData.mat');

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

sp1.PlotSpring(vecnorm([Qx(n), Qy(n)]), 'Configuration', [0, 0, atan2(Qy(n), Qx(n))]);

pause(.5);

t=0; %time initialization
n = 1;
tic;    %start time measuring
while (n <= numel(time))
    sp1.PlotSpring(vecnorm([Qx(n), Qy(n)]), 'Configuration', [0, 0, atan2(Qy(n), Qx(n))]);
    set(box_patch, 'XData', [0, 1.25, 1.25, 0] + Qx(n));
    drawnow;
    t = toc; %measuring current time
    n = round(t*freq)+1;
end
