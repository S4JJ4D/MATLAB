% RUN ME!

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
xline(0, 'k--');
yline(0, 'k--');
axis([-5 5 -5 5]);
ax.DataAspectRatio = [1 1 1];

%% Run Simulation
% loads pre-computed trajectory variables: Qx, Qy, time
load('motionData.mat');
freq = 1/median(diff(time));

bufferSize = 50;
ball_tl = TracingLine('ball_tl', bufferSize, ...
    'xBufferInit', Qx(1)*ones(1,bufferSize), ...
    'yBufferInit', Qy(1)*ones(1,bufferSize), ...
    'axes', ax);
ball_tl.Plot();

pause(1);

t=0; %time initialization
n = 1;
tic;    %start time measuring
while (n <= numel(time))

    ball_tl.AddPoint(Qx(n), Qy(n));
    drawnow;  %updates the display
    t = toc; %measuring current time
    n = round(t*freq)+1;
end
