function objects_vec = RendezvousCircle(ax, Ri, Ro, N, circle_center, options)
% RendezvousCircleBase  Draws Rendezvous Circle
%   objects_vec = RendezvousCircleBase(ax, Ri, Ro, N, circle_center,
%   options) draws a circle on the given axes with inner radius Ri, outer
%   radius Ro, divided into N sections and transported into the point
%   circle_center
%
% example: objects_vec = RendezvousCircleBase(ax, 1, 5, 8, [1, 2], 'scale', 2, 'arrows_count', 3);
arguments
    %
    ax (1,1)  {mustBeA(ax,["matlab.graphics.axis.Axes"," 'matlab.ui.control.UIAxes'"]), mustBeNonDeletedGraphicalObject}
    Ri              (1,1)  {mustBeNumeric,mustBeReal}
    Ro              (1,1)  {mustBeNumeric,mustBeReal}
    N               (1,1)  {mustBeNumeric,mustBeReal}
    circle_center   (1,2)  {mustBeNumeric,mustBeReal}
    options.outer_circle_color     (1,1)  string = "#502d16";
    options.inner_circle_color     (1,1)  string = "#aa8800";
    options.arrows_color     (1,1)  string = "#a02c2c";
    options.arrows_count     (1,1) {mustBePositive, mustBeInteger} = 1;
    options.scale (1,1) {mustBeNumeric,mustBeReal} = 1;
end


%%
% Draw Exterior Circle

c0 = circle_center;
outer_circle_color = options.outer_circle_color;
inner_circle_color = options.inner_circle_color;
arrows_color = options.arrows_color;
arrows_count = options.arrows_count;
scale_val = options.scale;

% 3+arrows_count gobjects are generated in this function
objects_vec = gobjects(1, 3+arrows_count);

hg = hgtransform(ax);
hg.Tag = 'RendezvousCircle_hg';
objects_vec(1) = hg;


t = 0:.01:2*pi;
co = plot(ax, Ro*cos(t), Ro*sin(t), 'LineWidth', 1.5, 'Color', outer_circle_color, ...
    'Parent', hg, ...
    'Tag', 'RendezvousCircle_OuterCircle');
objects_vec(2) =  co;

% Draw Interior Circle
ci = plot(ax, Ri*cos(t), Ri*sin(t), 'LineWidth', 1.5, 'Color', inner_circle_color, ...
    'Parent', hg, ...
    'Tag', 'RendezvousCircle_InnerCircle');
objects_vec(3) =  ci;

% divide the circle for arrows:

% sections for time
t_s = 0:2*pi/N:(2*pi - 2*pi/N);
discrete_points_out = [(Ro*cos(t_s));(Ro*sin(t_s))];
discrete_points_in = [(Ri*cos(t_s));(Ri*sin(t_s))];
vector_toward_center = discrete_points_in - discrete_points_out;

coeffs = (1/arrows_count):(1/arrows_count):1;

for i=1:arrows_count

    objects_vec(3+i) = ...
        quiver(ax, discrete_points_out(1,:), discrete_points_out(2,:), ...
        coeffs(i)*vector_toward_center(1,:), coeffs(i)*vector_toward_center(2,:), ...
        'AutoScale', 'off', 'LineWidth', 1.5, 'Color', arrows_color, ...
        'Parent', hg, ...
        'Tag', ['RendezvousCircle_Arrows_', num2str(i)]);

end

% qp_short = quiver(ax, discrete_points_out(1,:), discrete_points_out(2,:), vector_toward_center(1,:)/2, vector_toward_center(2,:)/2, ...
%     'AutoScale', 'off', 'LineWidth', 1.5, 'Color', arrows_color, ...
%     'Parent', hg, ...
%     'Tag', 'RendezvousCircle_ShortArrows');
% 
% qp_long = quiver(ax, discrete_points_out(1,:), discrete_points_out(2,:), vector_toward_center(1,:), vector_toward_center(2,:), ...
%     'AutoScale', 'off', 'LineWidth', 1.5, 'Color', arrows_color, ...
%     'Parent', hg, ...
%     'Tag', 'RendezvousCircle_LongArrows');

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