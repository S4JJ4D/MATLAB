%%
clear;
close all;
%%
r0 = .1; % roller radius
rb = .5; % base circle radius

l1 = 1.2; % link 1 length

bgColor = "#000000";
% cam1Color = "#212121";
cam1Color = "#1B222C";
linkColor = "#FFFFFF";
cam2Color = "#555555";

% cam2Color = '#9D7B53';

% baseCircleColor = "#BA0000";
baseCircleColor = "#F80000";



%%
% 7th order polynomial for rise/fall of the cam
Tmax = 10; % time length of a single motion cycle
TerminalVal = .5; % start/end value in a given motion cycle

time_vec = 0:.01:Tmax;
% specify the desired distance of the roller center from the cam pinning
% point
Profile1 = ProfileBuilder([0 TerminalVal 2 1.3 4 .6 8 1.1 12 .4 Tmax TerminalVal]);
% Profile2 = ProfileBuilder([0 TerminalVal 1 .8 2 .5 3 .8 4 .5 5 .8 6 .5 7 .8 8 .5 9 .8 Tmax TerminalVal]);
% Profile2 = ProfileBuilder([0 TerminalVal 3 1.2 4 .9 8 1.1 Tmax TerminalVal]);
Profile2 = ProfileBuilder([0 TerminalVal 6 1.6 Tmax TerminalVal]);
% Profile2 = ProfileBuilder([0 TerminalVal Tmax TerminalVal]);
[R_vec, R_dot_vec, R_ddot_vec] = Profile2(time_vec);

% Dependency of theta with time:
omega = 2*pi/Tmax;
h = @(x) omega*x;


fg_profile = figure;
set(fg_profile, 'Position', [633.8000 448.2000 560.0000 420], 'Units', 'pixels');

subplot(2,1,1);
plot(time_vec, R_vec, 'LineWidth', .8);
hold on;
grid on;
vline1 = xline(0, 'k-');
title('Motion Profile'); xlabel('time (s)'); ylabel('position (rad)');

subplot(2,1,2);
plot(time_vec, R_dot_vec, 'LineWidth', .8);
hold on;
vline2 = xline(0, 'k-');
grid on;
title('Time Derivative of Motion Profile'); xlabel('time (s)'); ylabel('velocity (rad/s)');

%%
% phi_vec1 = k(time_vec);
theta_vec = h(time_vec);


% Cam is rotating CCW
% E0 is the inertial frame located at the origin of the crack center of
% rotation

xc_E0 = R_vec;
yc_E0 = 0;

% xc is the relative position of the roller center w.r.t the cam
xc = R_vec.*cos(theta_vec);
yc = R_vec.*sin(theta_vec);

r_base = xc(1) - r0;

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
%         'HandleVisibility', 'on', 'Color', linkColor);
%     counter = counter + 1;
%     pause(.01);
% end

dR_dtheta = 1/omega * R_dot_vec;
d2R_dtheta2 = (1/omega)^2 * R_ddot_vec;

% Note that we parameterize the cam curve by theta.

% derivate with respect to theta
xc_p = dR_dtheta.*cos(theta_vec) - R_vec.*sin(theta_vec);
yc_p = dR_dtheta.*sin(theta_vec) + R_vec.*cos(theta_vec);

% second derivate with respect to theta
xc_pp = (d2R_dtheta2.*cos(theta_vec) - dR_dtheta.*sin(theta_vec)) - ...
        (dR_dtheta .*sin(theta_vec) + R_vec .* cos(theta_vec));
yc_pp = (d2R_dtheta2.*sin(theta_vec) + dR_dtheta.*cos(theta_vec)) + ...
        (dR_dtheta .*cos(theta_vec) - R_vec .* sin(theta_vec));



Delta = 1./((xc_p.^2 + yc_p.^2).^(1/2));
Delta_p = -(xc_p.^2 + yc_p.^2).^(-3/2) .* (xc_p .* xc_pp + yc_p .* yc_pp);

% xCam1 = xc + r0 * yc_p .* Delta;
% yCam1 = yc - r0 * xc_p .* Delta;

xCam2 = xc - r0 * yc_p .* Delta;
yCam2 = yc + r0 * xc_p .* Delta;

xCam2p = xc_p - r0 * (yc_pp .* Delta + yc_p .* Delta_p);
yCam2p = yc_p + r0 * (xc_pp .* Delta + xc_p .* Delta_p);

% now, that the cam x,y coordinates and its derivates are know, we can
% compute instantenous speed of the cam's curve
instantenous_speed_vec = sqrt(xCam2p.^2 + yCam2p.^2);
S = cumtrapz(theta_vec, instantenous_speed_vec);

% figure;
% subplot(2,1,1);
% plot(time_vec, instantenous_speed_vec);
% subplot(2,1,2);
% plot(time_vec, S);


%%
% cam1Plt = plot(xCam1, yCam1, 'b-', 'DisplayName', 'Cam1 Profile', 'LineWidth', 1.1);
% cam2Plt = plot(xCam2, yCam2, 'r-', 'DisplayName', 'Cam2 Profile', 'LineWidth', 1.1);
figure(fg_animation);
pause(.5);
cam2Plt = patch('XData', xCam2, 'YData', yCam2);
set(cam2Plt, 'FaceColor', cam2Color, 'EdgeColor', cam2Color, 'LineStyle', 'none');


%
% hold on;
% quiver(xCam2, yCam2, xCam2p, yCam2p, "AutoScaleFactor", .5);


%

% draw base circle with line
base_circle_plt = plot(r_base*cos(theta_vec),r_base*sin(theta_vec), 'Color', linkColor, 'DisplayName', 'Base Circle');
base_circle_line_plt = plot([0, xCam2(1)], [0, yCam2(1)], 'Color', linkColor);

% draw links
link1 = plot([0, l1], [0, 0], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', 2, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

% draw roller
roller_plt(1) = patch('XData', r0*cos(theta_vec), 'YData', r0*sin(theta_vec), ...
    'FaceColor', linkColor, 'FaceAlpha', 1);

numLines = 2;
for i=1:numLines
    roller_plt(1+i) = plot([0, r0*cos(i * 2*pi/numLines)], [0 r0*sin(i * 2*pi/numLines)], 'Color', 'k');
end

% roller_plt(2) = plot([0, -r0], [0 0], 'Color', 'k');


link1HTrans = hgtransform;
camHTrans = hgtransform;
rollerHTrans = hgtransform;

set(link1, 'Parent', link1HTrans, 'Color', linkColor, 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor);
set([cam2Plt, base_circle_plt, base_circle_line_plt] , 'Parent', camHTrans);
set(roller_plt, 'Parent', rollerHTrans);

% draw Cam1 Center
plot(0, 0, 'o', 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor, 'MarkerSize', 8);

link1HTrans.Matrix = makehgtform('translate', [R_vec(1), 0, 0]);
camHTrans.Matrix = makehgtform('zrotate', theta_vec(1));
rollerHTrans.Matrix = makehgtform('translate', [R_vec(1), 0, 0]) * makehgtform('zrotate', -S(1)/r0);

axis(ax, [-1.9550    2.9715   -1.6730    1.8723]);
ax.DataAspectRatio = [1 1 1];


%%
pause(1);
freq = 1/median(diff(time_vec));
t=0; %time initialization
n = 1;
tic;    %start time measuring
% I added the "ishandle" so the program will end in case u closed the figure
while (n <= numel(time_vec))
    link1HTrans.Matrix = makehgtform('translate', [R_vec(n), 0, 0]);
    camHTrans.Matrix = makehgtform('zrotate', -theta_vec(n));
    rollerHTrans.Matrix = makehgtform('translate', [R_vec(n), 0, 0]) * makehgtform('zrotate', S(n)/r0);

    set([vline1, vline2], 'Value', time_vec(n));
    drawnow;
    t = toc; %measuring current time
    n = round(t*freq)+1;
end