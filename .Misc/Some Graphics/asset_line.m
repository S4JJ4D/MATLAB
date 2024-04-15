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

% origin = plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

ps = [2, 1];
pe = [5, 8];

d = pe - ps;

points_num = 10;
coeffs = linspace(0, 1, points_num);

points = ps' + d' .* coeffs;
p = plot(points(1,:), points(2,:), 'Color', [0, 0, 1, 0.3], 'LineWidth', 1.5, 'Marker', '^');


% create patch objects





axis equal;

% y = linspace(0, 1, 10);
% x = 0 * ones(1,10);
% p = plot(x, y, 'Color', 'b', 'LineWidth', 1.5, 'Marker', '^');
% 
% hg = hgtransform;
% p.Parent = hg;
% 
% hg.Matrix = makehgtform('zrotate',pi/6);


