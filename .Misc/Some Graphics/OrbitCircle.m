function objects_vec = OrbitCircle(ax, Ri, Ro, circle_count, circle_center, options)
% OrbitCircle  Draws Orbiting Circles
%   Good scale value = 1/75 (: for Delta (x-axis length) = 25, use scale 1)
%
% objects_vec = OrbitCircle(ax, 1, 5, 4, [12, 20], 'scale', 2, 'arrow_count', 4)
arguments
    %
    ax (1,1)  {mustBeA(ax,["matlab.graphics.axis.Axes"," 'matlab.ui.control.UIAxes'"]), mustBeNonDeletedGraphicalObject}
    Ri              (1,1)  {mustBeNumeric,mustBeReal}
    Ro              (1,1)  {mustBeNumeric,mustBeReal}
    circle_count               (1,1)  {mustBePositive, mustBeInteger}
    circle_center   (1,2)  {mustBeNumeric,mustBeReal}
    options.circle_color     (1,1)  string = "#502d16";
    options.arrow_color (1,1) string = "#a02c2c";
    options.arrow_count     (1,1) {mustBeMember(options.arrow_count, [1, 2, 4])} = 1;
    options.scale (1,1) {mustBeNumeric,mustBeReal} = 1;
end


%%
% Draw Exterior Circle

c0 = circle_center;
circle_color = options.circle_color;
arrow_color = options.arrow_color;

arrow_count = options.arrow_count;
scale_val = options.scale;

% 4 gobjects are generated in this function
objects_vec = gobjects(1, circle_count+arrow_count);

hg = hgtransform(ax);
hg.Tag = 'OrbitCircle_hg';
j = 1;
objects_vec(j) = hg;



t = 0:.01:2*pi;
mrksize = 5;


radius_vec = linspace(Ri, Ro, circle_count);

circle_vec = gobjects(1, circle_count);
arrow_vec = gobjects(1, arrow_count);

for i=1:circle_count
    circle_vec(i) = plot(ax, radius_vec(i)*cos(t), radius_vec(i)*sin(t), ...
        'LineWidth', 1.5, 'Color', circle_color, 'Parent', hg, ...
        'Tag', ['OrbitCircle_Circle_', num2str(i)]);
    objects_vec(j+i) = circle_vec(i);
end
j = j + circle_count;

switch arrow_count
    case 1
        for i=1:arrow_count
            arrow_vec(i) = plot(ax, radius_vec*cos(pi/2), radius_vec*sin(pi/2), ...
                'Marker', '<', 'MarkerFaceColor', arrow_color, 'MarkerEdgeColor', arrow_color, ...
                'MarkerSize', mrksize, 'LineStyle', 'none', 'Parent', hg, ...
                'Tag', ['OrbitCircle_Arrow_', num2str(i)]);

            objects_vec(j+i) = arrow_vec(i);
        end

    case 2
        theta_vec = [pi/2, 3*pi/2];
        marker_vec = ['<', '>'];
        for i=1:arrow_count
            arrow_vec(i) = plot(radius_vec*cos(theta_vec(i)), radius_vec*sin(theta_vec(i)), ...
                'Marker', marker_vec(i), 'MarkerFaceColor', arrow_color, 'MarkerEdgeColor', arrow_color, ...
                'MarkerSize', mrksize, 'LineStyle', 'none', 'Parent', hg, ...
                'Tag', ['OrbitCircle_Arrow_', num2str(i)]);

            objects_vec(j+i) = arrow_vec(i);
        end

    case 4
        theta_vec = [0, pi/2, pi, 3*pi/2];
        marker_vec = ['^', '<', 'v', '>'];
        for i=1:arrow_count
            arrow_vec(i) = plot(radius_vec*cos(theta_vec(i)), radius_vec*sin(theta_vec(i)), ...
                'Marker', marker_vec(i), 'MarkerFaceColor', arrow_color, 'MarkerEdgeColor', arrow_color, ...
                'MarkerSize', mrksize, 'LineStyle', 'none', 'Parent', hg, ...
                'Tag', ['OrbitCircle_Arrow_', num2str(i)]);

            objects_vec(j+i) = arrow_vec(i);
        end

end


hg.Matrix = makehgtform('translate', [c0, 0]) * makehgtform('scale', scale_val);



end %// function

%% Custom validation function
function mustBeNonDeletedGraphicalObject(a)
if ~isgraphics(a)
    eidType = 'mustBeNonDeletedGraphicalObject:deletedObject';
    msgType = 'graphical input must be non-deleted';
    throwAsCaller(MException(eidType,msgType))
end
end

% 
% if exist('fg', 'var')
%     if ishandle(fg)
%         delete(fg);
%     end
% end
% clear;
% 
% Ro = 5;     % Outer Radius
% Ri = 1;     % Inner Radius
% N = 5;      % Resolution
% 
% c0 = [1, 2]; % circle center
% 
% 
% fg = figure(...
%     'Units', 'normalized', ...
%     'Position', [0.6453 0.3417 0.3464 0.5444]);
% ax = axes;
% hold(ax, 'on');
% grid on;
% box on;
% 
% xline(0, 'r--');
% yline(0, 'r--');
% 
% %%
% 
% t = 0:.01:2*pi;
% 
% hg = hgtransform; % used for scaling;
% 
% radius_vec = linspace(Ri, Ro, N);
% circle_vec = gobjects(1, N);
% arrow_vec = gobjects(1, N);
% for i=1:N
%     circle_vec(i) = plot(radius_vec(i)*cos(t) + c0(1), radius_vec(i)*sin(t) + c0(2), ...
%         'LineWidth', 1.5, 'Color', '#502d16', 'Parent', hg);
%     arrow_vec(i) = plot(radius_vec(i)*cos(pi/2) + c0(1), radius_vec(i)*sin(pi/2) + c0(2), ...
%         'Marker', '<', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
%         'MarkerSize', 7, 'Parent', hg);
% 
% end
% 
% 
% % co = plot(Ro*cos(t) + c0(1), Ro*sin(t) + c0(2), 'LineWidth', 1.5, 'Color', '#502d16');
% %
% % % Draw Interior Circle
% % ci = plot(Ri*cos(t) + c0(1), Ri*sin(t) + c0(2), 'LineWidth', 1.5, 'Color', '#aa8800');
% %
% % % divide the circle for arrows:
% %
% % % sections for time
% % t_s = 0:2*pi/N:(2*pi - 2*pi/N);
% % po = [(Ro*cos(t_s)+c0(1));(Ro*sin(t_s) + c0(2))];
% % pi = [(Ri*cos(t_s)+c0(1));(Ri*sin(t_s) + c0(2))];
% % vector_toward_center = pi - po;
% %
% % quiver(po(1,:), po(2,:), vector_toward_center(1,:), vector_toward_center(2,:), ...
% %     'AutoScale', 'off', 'LineWidth', 1.5, 'Color', '#a02c2c');
% 
% axis equal;

