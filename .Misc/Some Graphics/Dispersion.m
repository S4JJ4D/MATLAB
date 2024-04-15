
if exist('fg', 'var')
    if ishandle(fg)
        delete(fg);
    end
end
clear;

Ro = 5;     % Outer Radius
Ri = 1;     % Inner Radius
N = 8;      % Resolution

c0 = [1, 2]; % circle center


fg = figure(...
    'Units', 'normalized', ...
    'Position', [0.6453 0.3417 0.3464 0.5444]);
ax = axes;
hold(ax, 'on');
grid on;
box on;

xline(0, 'r--');
yline(0, 'r--');

%%
% Draw Exterior Circle

t = 0:.01:2*pi;
co = plot(Ro*cos(t) + c0(1), Ro*sin(t) + c0(2), 'LineWidth', 1.5, 'Color', '#502d16', ...
    'Visible', 'off');

% Draw Interior Circle
ci = plot(Ri*cos(t) + c0(1), Ri*sin(t) + c0(2), 'LineWidth', 1.5, 'Color', '#aa8800', ...
    'Visible', 'on');

% divide the circle for arrows:

% sections for time
t_s = 0:2*pi/N:(2*pi - 2*pi/N);
po = [(Ro*cos(t_s)+c0(1));(Ro*sin(t_s) + c0(2))];
pi = [(Ri*cos(t_s)+c0(1));(Ri*sin(t_s) + c0(2))];
vector_toward_center = po - pi;

quiver(pi(1,:), pi(2,:), vector_toward_center(1,:), vector_toward_center(2,:), ...
    'AutoScale', 'off', 'LineWidth', 1.5, 'Color', '#a02c2c');

axis equal;

