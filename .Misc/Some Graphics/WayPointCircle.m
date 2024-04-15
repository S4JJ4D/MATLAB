function objects_vec = WayPointCircle(ax, in_line, out_line, center, options)
% RendezvousCircleBase  Draws Rendezvous Circle
%   Good scale value = 1/25 (: for Delta (x-axis length) = 25, use scale 1)
%
% example: objects_vec4 = WayPointCircle(ax, [3, 1], [-2, 3], [0, 0], ...
%     'scale', 1, 'arrow_count', 3, ...
%     'offset', .3, 'arrow_color', '#4b184f');
arguments
    %
    ax (1,1)  {mustBeA(ax,["matlab.graphics.axis.Axes"," 'matlab.ui.control.UIAxes'"]), mustBeNonDeletedGraphicalObject}
    in_line              (1,2)  {mustBeNumeric,mustBeReal}
    out_line              (1,2)  {mustBeNumeric,mustBeReal}
    center               (1,2)  {mustBeNumeric,mustBeReal}

    options.line_length     (1,1)  {mustBeNumeric,mustBeReal} = 1;
    options.arrow_color     (1,1)  string = "#a02c2c";
    options.arc_color     (1,1)  string = "#aa8800"; %#502d16
    options.arrow_count     (1,1) {mustBePositive, mustBeInteger} = 1;
    options.dist_to_origin (1,1)  {mustBeNumeric,mustBeReal} = 1;
    options.offset (1,1)  {mustBeNumeric,mustBeReal} = 0;
    options.scale (1,1) {mustBeNumeric,mustBeReal} = 1;
end


%%
% Draw Exterior Circle

c0 = center;

line_length = options.line_length;
arrow_color = options.arrow_color;
arc_color = options.arc_color;
arrow_count = options.arrow_count;

dist_to_origin = options.dist_to_origin;
offset_val = options.offset;
scale_val = options.scale;

% 4 gobjects are generated in this function
objects_vec = gobjects(1, 2+2*arrow_count);

hg = hgtransform(ax);
hg.Tag = 'WayPointCircle_hg';
objects_vec(1) = hg;


% incoming vector
a = in_line;
% outgoing vector
b = out_line;

d1 = -a;
d2 = b;
d1_hat = d1/norm(d1);
d2_hat = d2/norm(d2);
e = d1_hat + d2_hat;
e_hat = e/norm(e);
% offset vector
offv = offset_val * e_hat;

delta = acos( dot(d1, d2)/(norm(d1) * norm(d2)) ) / 2;

t = 0:.01:2*pi;
% given desired radius
R = 1;
L = dist_to_origin;

givenR = 0;
% circumscribed circle center
if givenR
    c = R/sin(delta) * e_hat;
    L = R/sin(delta)
else
%     plot(L*cos(t), L*sin(t), 'k--');
    c = L * e_hat;
    R = L*sin(delta);
end



in_p2 = L*cos(delta) * d1_hat;
in_p1 = in_p2 + line_length * d1_hat;

in_dir = in_p2 - in_p1;

coeffs = (1/arrow_count):(1/arrow_count):1;

for i=1:arrow_count
objects_vec(1+i) = quiver(ax, offv(1)+in_p1(1), offv(2)+in_p1(2), ...
    coeffs(i)*in_dir(1), coeffs(i)*in_dir(2), 'AutoScale','off', ...
    'MaxHeadSize', .5, 'Color', arrow_color, ...
    'LineWidth', 2, ...
    'Parent', hg, ...
    'Tag', ['WayPointCircle_Arrow_in_', num2str(i)]);

end

out_p1 = L*cos(delta) * d2_hat;
out_p2 = out_p1 + line_length * d2_hat;

out_dir = out_p2 - out_p1;

for i=1:arrow_count
objects_vec(1+arrow_count+i) = quiver(ax, offv(1)+out_p1(1), offv(2)+out_p1(2), ...
    coeffs(i)*out_dir(1), coeffs(i)*out_dir(2), 'AutoScale','off', ...
    'MaxHeadSize', .5, 'Color', arrow_color, ...
    'LineWidth', 2, ...
    'Parent', hg, ...
    'Tag', ['WayPointCircle_Arrow_out', num2str(i)]);

end

r1 = in_p2 - c;
r2 = out_p1 - c;

theta_1 = atan2(r1(2), r1(1));
theta_2 = atan2(r2(2), r2(1));

t1 = min([theta_1, theta_2]);
t2 = max([theta_1, theta_2]);

if t2 - t1 > pi
    tt1 = wrapToPi(t2);
    tt2 = wrapToPi(t1);

else
    tt1 = t1;
    tt2 = t2;
end


t_sec = tt1:.01:tt2;

arc = plot(ax, R*cos(t_sec) + c(1) + offv(1), R*sin(t_sec) + c(2) + offv(2), ...
    'LineWidth', 2, 'Color', arc_color, ...
    'Parent', hg, ...
    'Tag', 'WayPointCircle_Arc');

objects_vec(end) = arc;

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