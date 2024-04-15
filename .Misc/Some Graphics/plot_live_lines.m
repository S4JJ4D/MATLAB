function [] = plot_live_lines(lines_list)
%PLOT_PROGRESSIVE_LIVE Summary of this function goes here
%   Detailed explanation goes here

% examples:
% lines_list = [[2;1], [5;8], [11;7], [12; 0]];
% plot_live_lines(lines_list)

% parameters:
% line_alpha: line alpha: 
% points_num: number of points per line
% arrows_num: number of arrows per line



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

lines_num = size(lines_list, 2) - 1;

data_struct = repmat(struct(), 1, lines_num);

line_alpha = .2;
p_whole_line = plot(lines_list(1,:), lines_list(2,:), 'Color', [1, 0, 0, line_alpha], ...
        'LineWidth', 6, 'Marker', '.', 'MarkerSize', 4, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b');
for i=1:lines_num

    d = lines_list(:,i+1) - lines_list(:,i);
    points_num = 17;
    coeffs = linspace(0, 1, points_num);

    points = lines_list(:,i) + d .* coeffs;
    p = plot(points(1,:), points(2,:), 'Color', [1, 0, 0, 0.0], ...
        'LineWidth', 6, 'Marker', '.', 'MarkerSize', 4, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b');
    line_ang = atan2(d(2), d(1));

    data_struct(i).points = points;
    data_struct(i).line_plt = p;
    data_struct(i).line_ang = line_ang;

end


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

for j=1:lines_num

    points = data_struct(j).points;
    line_ang = data_struct(j).line_ang;

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

    data_struct(j).arrows_list = arrows_list;
    data_struct(j).hg_list = hg_list;

end
% axis equal;
axis([-3.6114   18.0198   -5.0875   15.0219]);
% animate

% frequency
w = .8;
time = 1:.1:20;
% phase lage
alpha = linspace(0,1.5*pi, arrows_num+10);
% alpha = linspace(0, 2*pi, lines_num*arrows_num);

for t=time
    for j=1:lines_num
        for i=1:arrows_num
            data_struct(j).arrows_list(i).FaceAlpha = ...
                1/2 + 1/2*sin(w*(t - alpha(i)));
%              1/2 + 1/2*sin(w*(t - alpha(i+(j-1)*arrows_num)));
        end
        drawnow;
    end
    pause(.02);
end




end

