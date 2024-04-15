% RUN ME!

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

points_num = 17;
coeffs = linspace(0, 1, points_num);

points = ps' + d' .* coeffs;
p = plot(points(1,:), points(2,:), 'Color', [0, 0, 1, 0.2], ...
    'LineWidth', 5, 'Marker', '.', 'MarkerSize', 8);

line_ang = atan2(d(2), d(1));
% create patch objects
init_scl = .25;
x = 1.5;
y = 1/1.5 * x;
w = 1/2 * x;

X = [x, x, 0, -x, -x, 0];
Y = [0, w, y+w, w, 0, y];


% choose an odd number to preserve symmetry
arrows_num = 11;
if ~mod(arrows_num,2)
    arrows_num = arrows_num + 1;
end
k = (arrows_num - 1)/2;

idx_range = (ceil(points_num/2) - k):(ceil(points_num/2) + k);

arrows_list = gobjects(1, arrows_num);

hg_list =  gobjects(1, arrows_num);
for i=1:arrows_num
    hg_list(i) = hgtransform;
    arrows_list(i) = patch('XData',X,'YData',Y, ...
        'Parent', hg_list(i), 'LineStyle', 'none');
    hg_list(i).Matrix = ...
        makehgtform('translate', [points(1, idx_range(i)), points(2, idx_range(i)), 0]) * ...
        makehgtform('scale',init_scl) * ...
        makehgtform('zrotate', -(pi/2 - line_ang));
end

axis equal;

% animate

% frequency
w = .7;
time = 1:.1:20;
% phase lage
alpha = linspace(0, 1.5*pi, arrows_num+5);
% alpha = linspace(0, 2*pi, arrows_num);

for t=time
    for i=1:arrows_num
        arrows_list(i).FaceAlpha = 1/2 + 1/2*sin(w*(t - alpha(i)));
    end
    drawnow;
    pause(.02);
end



