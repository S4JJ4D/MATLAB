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

x = 1.5;
y = 1/1.5 * x;
w = 1/2 * x;

X = [x, x, 0, -x, -x, 0];
Y = [0, w, y+w, w, 0, y];

arr_p = patch('XData',X,'YData',Y);


% origin = plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
% create patch objects

hg = hgtransform;
arr_p.Parent = hg;

hg.Matrix = makehgtform('zrotate',pi/6);


axis equal;

% y = linspace(0, 1, 10);
% x = 0 * ones(1,10);
% p = plot(x, y, 'Color', 'b', 'LineWidth', 1.5, 'Marker', '^');
% 



