clear all;
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
xline(0, 'r--');
yline(0, 'r--');
axis([ -4.7956    4.8454   -4.9919    5.6491]);
ax.DataAspectRatio = [1 1 1];

sp1 = DynamicSpring('spring1', 'Radius', .25, 'Pitch', .03, 'Turns', 18, 'Axes', ax);
sp2 = DynamicSpring('spring2', 'Radius', .25, 'Pitch', .03, 'Turns', 13, 'Axes', ax);
sp2.plotting_options.FrontEyeOuterColor = 'b';
sp2.plotting_options.RearEyeOuterColor = 'r';


r_ball = .2;
theta_vec = 0:.1:2*pi;
ball1_plt = patch('XData', r_ball*cos(theta_vec), 'YData', r_ball*sin(theta_vec));
ball2_plt = patch('XData', r_ball*cos(theta_vec), 'YData', r_ball*sin(theta_vec));

hg_ball1 = hgtransform;
hg_ball2 = hgtransform;


set(hg_ball1, 'Tag', 'ball1');
set(hg_ball2, 'Tag', 'ball2');

set(ball1_plt, 'Parent', hg_ball1, 'FaceAlpha', 1, 'FaceColor', 'b');
set(ball2_plt, 'Parent', hg_ball2, 'FaceAlpha', 1, 'FaceColor', 'r');

% trace1_plt = plot(0, 0, 'b--');
% trace2_plt = plot(0, 0, 'r--');
% set([trace1_plt, trace2_plt], 'XData', [], 'YData', []);

%% Run Simulation
% simout = sim('DoubleMassSimulaitonModel');
load('double_mass_spring.mat');
% time = simout.tout.';
freq = 1/median(diff(time));

% q_vec = squeeze(simout.logsout.getElement('q').Values.Data);
% x1 = q_vec(1,:);
% y1 = q_vec(2,:);
% x2 = q_vec(3,:);
% y2 = q_vec(4,:);

n = 1;

% sp1.plotting_options.LineWidth = 2.5;
% sp2.plotting_options.LineWidth = 2.5;

sp1.PlotSpring(vecnorm([x1(n), y1(n)]), 'Configuration', [0, 0, atan2(y1(n), x1(n))]);
sp2.PlotSpring(vecnorm([x2(n), y2(n)]), 'Configuration', [x1(n), y1(n), atan2(y2(n), x2(n))]);


% sp1.PlotSpring([0, 0], atan2(y1(n), x1(n)), vecnorm([x1(n), y1(n)]), axes=ax, color='#0056A1');
% sp2.PlotSpring([x1(n), y1(n)], atan2(y2(n), x2(n)), vecnorm([x2(n), y2(n)]), axes=ax, color='#4a691b');

hg_ball1.Matrix = makehgtform('translate', [x1(n), y1(n), 0]);
hg_ball2.Matrix = makehgtform('translate', [x1(n)+x2(n), y1(n)+y2(n), 0]);

set([hg_ball1, hg_ball2], 'Visible', 'off');

% bring balls to the top
ax.Children = circshift(ax.Children, -2);

pause(1);

t=0; %time initialization
n = 1;
tic;    %start time measuring
% I added the "ishandle" so the program will end in case u closed the figure
while (n <= numel(time))

    %     sp1.PlotSpring([0, 0], atan2(y1(n), x1(n)), vecnorm([x1(n), y1(n)]));
    %     sp2.PlotSpring([x1(n), y1(n)], atan2(y2(n), x2(n)), vecnorm([x2(n), y2(n)]));

    sp1.PlotSpring(vecnorm([x1(n), y1(n)]), 'Configuration', [0, 0, atan2(y1(n), x1(n))]);
    sp2.PlotSpring(vecnorm([x2(n), y2(n)]), 'Configuration', [x1(n), y1(n), atan2(y2(n), x2(n))]);


    hg_ball1.Matrix = makehgtform('translate', [x1(n), y1(n), 0]);
    hg_ball2.Matrix = makehgtform('translate', [x1(n)+x2(n), y1(n)+y2(n), 0]);

    %     set(trace1_plt, 'XData', [trace1_plt.XData, x1(n)], 'YData', [trace1_plt.YData, y1(n)]);
    %     set(trace2_plt, 'XData', [trace2_plt.XData, x1(n)+x2(n)], 'YData', [trace2_plt.YData, y1(n)+y2(n)]);


    drawnow;  %updates the display
    t = toc; %measuring current time
    n = round(t*freq)+1;
end
