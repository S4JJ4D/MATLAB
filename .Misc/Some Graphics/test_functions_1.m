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

objects_vec = RendezvousCircle(ax, 1, 5, 8, [20, 10], 'scale', 2, 'arrows_count', 3);
objects_vec2 = DispersionCircle(ax, 1, 5, 8, [-25, 0], 'scale', 2, 'arrows_count', 3);
objects_vec3 = OrbitCircle(ax, 1, 5, 4, [12, -20], 'scale', 2, 'arrow_count', 2);
objects_vec4 = WayPointCircle(ax, [3, 1], [-2, 3], [4, 3], 'scale', 11, 'arrow_count', 3);



axis(ax, 'equal');