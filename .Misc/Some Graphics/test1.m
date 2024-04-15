% RUN ME!

if exist('fg', 'var')
    if ishandle(fg)
        delete(fg);
    end
end
clear;
fg = figure(...
    'Units', 'normalized', ...
    'Position', [0.6453 0.3417 0.3464 0.5444]);
ax = axes;
hold(ax, 'on');
grid on;
box on;

xline(0, 'r--');
yline(0, 'r--');

% objects_vec4 = WayPointCircle(ax, [3, 1], [-2, 3], [0, 0], ...
%     'scale', 1, 'arrow_count', 3, ...
%     'offset', .3, 'arrow_color', '#4b184f');

objects_vec3 = OrbitCircle(ax, 1, 5, 4, [0, 0], 'scale', 1, 'arrow_count', 1);




axis(ax, 'equal');