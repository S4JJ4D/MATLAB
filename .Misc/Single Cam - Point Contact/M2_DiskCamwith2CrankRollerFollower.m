%%
clear;
close all;
%%
rc = 0.8;   % although, each cam could have its own rc
r0 = .05; % although, each cam could have its own r0

l1 = 1.2; % link 1 length
a1 = .6;
b1 = .3;

l2 = 2.3;
a2 = .35;
b2 = .25;

markerSize = 5;

bgColor = "#000000";
% cam1Color = "#212121";
cam1Color = "#1B222C";
linkColor = "#FFFFFF";
% linkColor = "#000000";
cam2Color = "#555555";

% baseCircleColor = "#BA0000";
baseCircleColor = "#F80000";


drawStar = true;

%% Arbitrary Profile
% 7th order polynomial for rise/fall of the cam
if drawStar == false

    Tmax = 15; % time length of a single motion cycle
    TermVal1 = .5; % start/end value in a given motion cycle
    TermVal2 = -.6; % start/end value in a given motion cycle
    time_vec = 0:.01:Tmax;
    Profile1 = ProfileBuilder([0 TermVal1 6 1.3 12 1 Tmax TermVal1]);
    Profile2 = ProfileBuilder([0 TermVal2 3 -.8 Tmax TermVal2]);
    [phi_vec1, phi_dot_vec1] = Profile1(time_vec);
    [phi_vec2, phi_dot_vec2] = Profile2(time_vec);

else
    %% Linear Motion on the task-space with Inverse Kinematics
    Tmax = 10;
    time_vec = 0:.01:Tmax;
    pMat = [[2.5, -.5];[3, .1];[2.5, .3]];

    star_coord = [...
        [1,0.272727272727273];
        [0.238471673254282,0.272727272727273];
        [0.001317523056654,1];
        [-0.235836627140975,0.275362318840580];
        [-1,0.272727272727273];
        [-0.383399209486166,-0.175230566534914];
        [-0.617918313570488,-0.899868247694335];
        [0.003952569169960,-0.454545454545455];
        [0.617918313570488,-0.899868247694335];
        [0.383399209486166,-0.175230566534914]];

    ps = star_coord.';
    ScalingMat = .3*eye(2);
    ps2 = ScalingMat * ps;

    tt = 15;
    RotationMat = [cosd(tt), -sind(tt);sind(tt), cosd(tt)];

    ps3 = RotationMat * ps2 + [2.52; -0];
    pMat = ps3.';
    LineProfile = LinearPathProfile(pMat, l1, l2, Tmax);
    [phi_vec1, phi_dot_vec1, phi_vec2, phi_dot_vec2] = LineProfile(time_vec);

end
%%

% Dependency of theta with time:
omega = 2*pi/Tmax;
h = @(x) omega*x;


fg_profile1 = figure;
set(fg_profile1, 'Position', [630.6000 473 560.0000 402.4000], 'Units', 'pixels');

subplot(2,1,1);
plot(time_vec, phi_vec1, 'LineWidth', .8);
hold on;
grid on;
vline11 = xline(0, 'k-');
title('Motion Profile'); xlabel('time (s)'); ylabel('position (rad)');

subplot(2,1,2);
plot(time_vec, phi_dot_vec1, 'LineWidth', .8);
hold on;
vline12 = xline(0, 'k-');
grid on;
title('Time Derivative of Motion Profile'); xlabel('time (s)'); ylabel('velocity (rad/s)');


fg_profile2 = figure;
set(fg_profile2, 'Position', [631.4000 3.4000 560 388], 'Units', 'pixels');

subplot(2,1,1);
plot(time_vec, phi_vec2, 'LineWidth', .8);
hold on;
grid on;
vline21 = xline(0, 'k-');
title('Motion Profile'); xlabel('time (s)'); ylabel('position (rad)');

subplot(2,1,2);
plot(time_vec, phi_dot_vec2, 'LineWidth', .8);
hold on;
vline22 = xline(0, 'k-');
grid on;
title('Time Derivative of Motion Profile'); xlabel('time (s)'); ylabel('velocity (rad/s)');




%%
theta_vec1 = h(time_vec);

% Cam is rotating CCW
% E0 is the inertial frame located at the origin of the crack center of
% rotation

xc1_E0 = a1*cos(phi_vec1) + b1*sin(phi_vec1);
yc1_E0 = a1*sin(phi_vec1) - b1*cos(phi_vec1);

% xc is the relative position of the roller center w.r.t the cam
xc1 = a1*cos(phi_vec1 - theta_vec1) + b1*sin(phi_vec1 - theta_vec1) - rc*cos(theta_vec1);
yc1 = a1*sin(phi_vec1 - theta_vec1) - b1*cos(phi_vec1 - theta_vec1) + rc*sin(theta_vec1);

r_base1 = vecnorm([xc1(1), yc1(1)]) - r0;

fg_animation = figure;
set(fg_animation, 'Position', [67.4000 473.8000 560 401.6000], 'Units', 'pixels');
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

dphi_dtheta1 = 1/omega * phi_dot_vec1;

xc1_p = (-a1*sin(phi_vec1 - theta_vec1) + b1*cos(phi_vec1 - theta_vec1)).* (dphi_dtheta1 - 1) + rc*sin(theta_vec1);
yc1_p = (a1*cos(phi_vec1 - theta_vec1) + b1*sin(phi_vec1 - theta_vec1)) .* (dphi_dtheta1 - 1) + rc*cos(theta_vec1);

Delta1 =1./((xc1_p.^2 + yc1_p.^2).^(1/2));

xCam1 = xc1 + r0 * yc1_p .* Delta1;
yCam1 = yc1 - r0 * xc1_p .* Delta1;

%%
theta_vec2 = h(time_vec);


% Cam is rotating CCW
% E0 is the inertial frame located at the origin of the crack center of
% rotation

xc2_E0 = l1*cos(phi_vec1) + a2*cos(phi_vec1 + phi_vec2) + b2*sin(phi_vec1 + phi_vec2);
yc2_E0 = l1*sin(phi_vec1) + a2*sin(phi_vec1 + phi_vec2) - b2*cos(phi_vec1 + phi_vec2);

% xc is the relative position of the roller center w.r.t the cam
xc2 = a2*cos(phi_vec1 + phi_vec2 - theta_vec2) + b2*sin(phi_vec1 + phi_vec2 - theta_vec2) + ...
    -rc*cos(theta_vec2) + l1*cos(phi_vec1 - theta_vec2);

yc2 = a2*sin(phi_vec1 + phi_vec2 - theta_vec2) - b2*cos(phi_vec1 + phi_vec2 - theta_vec2) + ...
    +rc*sin(theta_vec2) + l1*sin(phi_vec1 - theta_vec2);

r_base2 = vecnorm([xc2(1), yc2(1)]) - r0;

% fg_animation = figure;
% set(fg_animation, 'Position', [67.4000 473.8000 560 401.6000], 'Units', 'pixels');
% plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
% ax = gca;
% % ax.Color = bgColor;
% hold on;
%
% axis equal;
% pause(1);
%
% i=1:10:numel(time_vec);
% rollers_center_plt_array = gobjects(1, numel(i));
% counter = 1;
% for i=1:10:numel(time_vec)
%     rollers_center_plt_array(counter) = plot(xc2(i) + r0*cos(theta_vec2), yc2(i) + r0*sin(theta_vec2), ...
%         'HandleVisibility', 'on', 'Color', 'k');
%     counter = counter + 1;
%     pause(.01);
% end

dphi_dtheta2 = 1/omega * phi_dot_vec2;

xc2_p = (-a2*sin(phi_vec1 + phi_vec2 - theta_vec2) + b2*cos(phi_vec1 + phi_vec2 - theta_vec2)) .* ...
    (dphi_dtheta1 + dphi_dtheta2 - 1) + ...
    rc*sin(theta_vec2) - l1*sin(phi_vec1 - theta_vec2) .* (dphi_dtheta1 - 1);

yc2_p = (a2*cos(phi_vec1 + phi_vec2 - theta_vec2) + b2*sin(phi_vec1 + phi_vec2 - theta_vec2)) .* ...
    (dphi_dtheta1 + dphi_dtheta2 - 1) + ...
    rc*cos(theta_vec2) + l1*cos(phi_vec1 - theta_vec2) .* (dphi_dtheta1 - 1);

Delta2 =1./((xc2_p.^2 + yc2_p.^2).^(1/2));

xCam2 = xc2 + r0 * yc2_p .* Delta2;
yCam2 = yc2 - r0 * xc2_p .* Delta2;

% xCam22 = xc2 - r0 * yc2_p .* Delta2;
% yCam22 = yc2 + r0 * xc2_p .* Delta2;

% cam1Plt = plot(xCam22, yCam22, 'b-', 'DisplayName', 'Cam1 Profile', 'LineWidth', 1.1);
% cam2Plt = plot(xCam2, yCam2, 'r-', 'DisplayName', 'Cam2 Profile', 'LineWidth', 1.1);


%%
% cam1Plt = plot(xCam1, yCam1, 'b-', 'DisplayName', 'Cam1 Profile', 'LineWidth', 1.1);
% cam2Plt = plot(xCam2, yCam2, 'r-', 'DisplayName', 'Cam2 Profile', 'LineWidth', 1.1);

cam2Plt = patch(ax, 'XData', xCam2, 'YData', yCam2, 'EdgeColor', 'b', 'LineWidth', .5, 'FaceAlpha', .2);
hold on;
cam1Plt = patch(ax, 'XData', xCam1, 'YData', yCam1, 'EdgeColor', 'r', 'LineWidth', .5, 'FaceAlpha', .2);

set(cam2Plt, 'FaceColor', cam2Color, 'EdgeColor', cam2Color, 'FaceAlpha', 1);
set(cam1Plt, 'FaceColor', cam1Color, 'EdgeColor', cam1Color, 'FaceAlpha', 1);

% draw base circle with line
% base_circle_plt = plot(r_base1*cos(theta_vec1),r_base1*sin(theta_vec1), 'Color', baseCircleColor, 'DisplayName', 'Base Circle');
% base_circle_line_plt = plot([0, xCam1(1)], [0, yCam1(1)], 'Color', baseCircleColor);

% draw links
link1(1) = plot([0, l1], [0, 0], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [1, 2], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', markerSize);
link1(2) = plot([a1, a1], [0, -b1], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', markerSize);

link2(1) = plot([0, l2], [0, 0], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', 1, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', markerSize);
link2(2) = plot([a2, a2], [0, -b2], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', markerSize);


% plot star:
if drawStar
    patch('XData', pMat(:,1), 'YData', pMat(:,2), 'FaceColor', 'g', 'FaceAlpha', '0', ...
        'EdgeColor', '#FCE500', 'LineWidth', 1.5);
end
% plot(pMat(:,1), pMat(:,2), 'y-');


cam2HTrans = hgtransform;
cam1HTrans = hgtransform;
link2HTrans = hgtransform;
link1HTrans = hgtransform;


set(link1, 'Parent', link1HTrans, 'Color', linkColor, 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor);
set(link2, 'Parent', link2HTrans, 'Color', linkColor, 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor);

% set([cam1Plt, base_circle_plt, base_circle_line_plt] , 'Parent', cam1HTrans);
set(cam2Plt , 'Parent', cam2HTrans);
set(cam1Plt , 'Parent', cam1HTrans);


% draw Cam1 Center
plot(rc, 0, 'o', 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor, 'MarkerSize', 8);

link1HTrans.Matrix = makehgtform('zrotate', phi_vec1(1));
link2HTrans.Matrix = makehgtform('zrotate', phi_vec1(1)) * makehgtform('translate', [l1 0 0]) * makehgtform('zrotate', phi_vec2(1));
cam1HTrans.Matrix = makehgtform('translate', [rc 0 0]) * makehgtform('zrotate', theta_vec1(1));
cam2HTrans.Matrix = makehgtform('translate', [rc 0 0]) * makehgtform('zrotate', theta_vec2(1));

% draw roller
roller1_plt = patch('XData', xc1_E0(1) + r0*cos(theta_vec1), 'YData', yc1_E0(1) + r0*sin(theta_vec1), ...
    'FaceColor', linkColor, 'FaceAlpha', 1);

roller2_plt = patch('XData', xc2_E0(1) + r0*cos(theta_vec2), 'YData', yc2_E0(1) + r0*sin(theta_vec2), ...
    'FaceColor', linkColor, 'FaceAlpha', 1);

axis(ax, [-0.2713    3.8150   -1.4480    1.7750]);
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
    link2HTrans.Matrix = makehgtform('zrotate', phi_vec1(n)) * makehgtform('translate', [l1 0 0]) * makehgtform('zrotate', phi_vec2(n));

    cam1HTrans.Matrix = makehgtform('translate', [rc 0 0]) * makehgtform('zrotate', theta_vec1(n));
    cam2HTrans.Matrix = makehgtform('translate', [rc 0 0]) * makehgtform('zrotate', theta_vec2(n));

    set([vline11, vline12, vline21, vline22], 'Value', time_vec(n));

    set(roller1_plt, 'XData', xc1_E0(n) + r0*cos(theta_vec1), ...
        'YData', yc1_E0(n) + r0*sin(theta_vec1));

    set(roller2_plt, 'XData', xc2_E0(n) + r0*cos(theta_vec2), ...
        'YData', yc2_E0(n) + r0*sin(theta_vec2));
    drawnow;
    t = toc; %measuring current time
    n = round(t*freq)+1;
end

