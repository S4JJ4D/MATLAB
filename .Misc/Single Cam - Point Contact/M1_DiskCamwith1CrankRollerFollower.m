%%
clear;
close all;
%%
rc = 1;
r0 = .05;

l1 = 1.2; % link 1 length
a1 = .8;
b1 = .25;


bgColor = "#000000";
% cam1Color = "#212121";
cam1Color = "#1B222C";
linkColor = "#FFFFFF";
cam2Color = "#555555";

% baseCircleColor = "#BA0000";
baseCircleColor = "#F80000";



%%
% 7th order polynomial for rise/fall of the cam
Tmax = 15; % time length of a single motion cycle
TermVal = .5; % start/end value in a given motion cycle

time_vec = 0:.01:Tmax;
Profile1 = ProfileBuilder([0 TermVal 2 1.3 4 .6 8 1.1 12 .4 Tmax TermVal]);
Profile2 = ProfileBuilder([0 TermVal 6 1.3 12 1 Tmax TermVal]);
[phi_vec1, phi_dot_vec1] = Profile2(time_vec);

% Dependency of theta with time:
omega = 2*pi/Tmax;
h = @(x) omega*x;


fg_profile = figure;
set(fg_profile, 'Position', [633.8000 448.2000 560.0000 420], 'Units', 'pixels');

subplot(2,1,1);
plot(time_vec, phi_vec1, 'LineWidth', .8);
hold on;
grid on;
vline1 = xline(0, 'k-');
title('Motion Profile'); xlabel('time (s)'); ylabel('position (rad)');

subplot(2,1,2);
plot(time_vec, phi_dot_vec1, 'LineWidth', .8);
hold on;
vline2 = xline(0, 'k-');
grid on;
title('Time Derivative of Motion Profile'); xlabel('time (s)'); ylabel('velocity (rad/s)');

%%
% phi_vec1 = k(time_vec);
theta_vec1 = h(time_vec);


% Cam is rotating CCW
% E0 is the inertial frame located at the origin of the crack center of
% rotation

xc_E0 = a1*cos(phi_vec1) + b1*sin(phi_vec1);
yc_E0 = a1*sin(phi_vec1) - b1*cos(phi_vec1);

% xc is the relative position of the roller center w.r.t the cam
xc = a1*cos(phi_vec1 - theta_vec1) + b1*sin(phi_vec1 - theta_vec1) - rc*cos(theta_vec1);
yc = a1*sin(phi_vec1 - theta_vec1) - b1*cos(phi_vec1 - theta_vec1) + rc*sin(theta_vec1);

r_base = vecnorm([xc(1), yc(1)]) - r0;

fg_animation = figure;
set(fg_animation, 'Position', [68.2000 449 560 420], 'Units', 'pixels');
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
ax = gca;
ax.Color = bgColor;
hold on;

axis equal;
pause(1);

% i=1:10:numel(time_vec);
% rollers_center_plt_array = gobjects(1, numel(i));
% counter = 1;
% for i=1:10:numel(time_vec)
%     rollers_center_plt_array(counter) = plot(xc(i) + r0*cos(theta_vec), yc(i) + r0*sin(theta_vec), ...
%         'HandleVisibility', 'on', 'Color', 'k');
%     counter = counter + 1;
%     pause(.01);
% end

dphi_dtheta = 1/omega * phi_dot_vec1;

xc_p = (-a1*sin(phi_vec1 - theta_vec1) + b1*cos(phi_vec1 - theta_vec1)).* (dphi_dtheta - 1) + rc*sin(theta_vec1);
yc_p = (a1*cos(phi_vec1 - theta_vec1) + b1*sin(phi_vec1 - theta_vec1)) .* (dphi_dtheta - 1) + rc*cos(theta_vec1);

Delta =1./((xc_p.^2 + yc_p.^2).^(1/2));

xCam1 = xc + r0 * yc_p .* Delta;
yCam1 = yc - r0 * xc_p .* Delta;

% xCam2 = xc - r0 * yc_p .* Delta;
% yCam2 = yc + r0 * xc_p .* Delta;

%%
% cam1Plt = plot(xCam1, yCam1, 'b-', 'DisplayName', 'Cam1 Profile', 'LineWidth', 1.1);
% cam2Plt = plot(xCam2, yCam2, 'r-', 'DisplayName', 'Cam2 Profile', 'LineWidth', 1.1);

cam1Plt = patch('XData', xCam1, 'YData', yCam1);
set(cam1Plt, 'FaceColor', cam2Color, 'EdgeColor', cam2Color);

% draw base circle with line
base_circle_plt = plot(r_base*cos(theta_vec1),r_base*sin(theta_vec1), 'Color', baseCircleColor, 'DisplayName', 'Base Circle');
base_circle_line_plt = plot([0, xCam1(1)], [0, yCam1(1)], 'Color', baseCircleColor);

% draw links
link1(1) = plot([0, l1], [0, 0], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [1, 2], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
link1(2) = plot([a1, a1], [0, -b1], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 7);

% draw roller
roller_plt = patch('XData', xc_E0(1) + r0*cos(theta_vec1), 'YData', yc_E0(1) + r0*sin(theta_vec1), ...
    'FaceColor', linkColor, 'FaceAlpha', 1);

link1HTrans = hgtransform;
camHTrans = hgtransform;
set(link1, 'Parent', link1HTrans, 'Color', linkColor, 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor);
set([cam1Plt, base_circle_plt, base_circle_line_plt] , 'Parent', camHTrans);

% draw Cam1 Center
plot(rc, 0, 'o', 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor, 'MarkerSize', 8);

link1HTrans.Matrix = makehgtform('zrotate', phi_vec1(1));
camHTrans.Matrix = makehgtform('translate', [rc 0 0]) * makehgtform('zrotate', theta_vec1(1));

axis(ax, [-0.8187, 2.8962, -1.4460, 1.4840]);
ax.DataAspectRatio = [1 1 1];


%%
pause(1);
freq = 1/median(diff(time_vec));
t=0; %time initialization
n = 1;
tic;    %start time measuring
% I added the "ishandle" so the program will end in case u closed the figure
while (n <= numel(time_vec))
    link1HTrans.Matrix = makehgtform('zrotate', phi_vec1(n));
    camHTrans.Matrix = makehgtform('translate', [rc 0 0]) * makehgtform('zrotate', theta_vec1(n));

    set([vline1, vline2], 'Value', time_vec(n));
    set(roller_plt, 'XData', xc_E0(n) + r0*cos(theta_vec1), ...
        'YData', yc_E0(n) + r0*sin(theta_vec1));
    drawnow;
    t = toc; %measuring current time
    n = round(t*freq)+1;
end

