%%
clear;
close all;
%%
% first component of the following vectors refer to the horizontal cam and
% the 2nd one refers to the vertical cam.
rf = [.12, .15]; % roller radius
rb = [.5, .8]; % base circle radius
e = [-.25, .4];  % offset

Rp = rf + rb; % called "prime circle radius"
beta = asin(e./(Rp)); % called "initial offset angle"
d = Rp.*cos(beta); % initial follower position at low dwel

l = [1.6, 1.2];
r_BA = [-.5, 2];



bgColor = "#000000";
% cam1Color = "#212121";
cam1Color = "#1B222C";
linkColor = "#FFFFFF";
cam2Color = "#555555";

vCamColor = '#705549';

% cam2Color = '#9D7B53';

% baseCircleColor = "#BA0000";
baseCircleColor = "#F80000";

%%
% 7th order polynomial for rise/fall of the cam
Tmax = 10; % time length of a single motion cycle
TerminalVal = 0; % start/end value in a given motion cycle

hh = Tmax/6;

time_vec = 0:.01:Tmax;
% specify the desired distance of the roller center from the cam pinning
% point
a = .5;
b = .25;
xProfile = ProfileBuilder([0 0 hh 0 2*hh 0 3*hh b 4*hh b 5*hh b 6*hh 0]);
yProfile = ProfileBuilder([0 0 hh a 2*hh a/2 3*hh a/2 4*hh 0 5*hh a 6*hh 0]);

[sx_vec, sx_dot_vec, sx_ddot_vec] = xProfile(time_vec);
[sy_vec, sy_dot_vec, sy_ddot_vec] = yProfile(time_vec);


% Dependency of theta with time:
omega = 2*pi/Tmax;
h = @(x) omega*x;


fg_profile = figure;
set(fg_profile, 'Position', [108 371 1124 583], 'Units', 'pixels', ...
    'Color', 'w');

tiledlayout(3,2,"TileSpacing",'compact','Padding','compact');
nexttile(2);
plot(time_vec, sx_vec, 'LineWidth', .8);
hold on;
grid on;
vline(1) = xline(0, 'k-');
title('Motion Profile: X'); xlabel('time (s)'); ylabel('position');

nexttile(4);
plot(time_vec, sx_dot_vec, 'LineWidth', .8);
hold on;
vline(2) = xline(0, 'k-');
grid on;
title('Time Derivative of Motion Profile: Xdot'); xlabel('time (s)'); ylabel('velocity (/s)');


nexttile(5);
plot(time_vec, sy_vec, 'LineWidth', .8);
hold on;
grid on;
vline(3) = xline(0, 'k-');
title('Motion Profile: Y'); xlabel('time (s)'); ylabel('position');

nexttile(6);
plot(time_vec, sy_dot_vec, 'LineWidth', .8);
hold on;
vline(4) = xline(0, 'k-');
grid on;
title('Time Derivative of Motion Profile: Ydot'); xlabel('time (s)'); ylabel('velocity (/s)');



%%
% phi_vec1 = k(time_vec);
theta_vec = h(time_vec);


% Cam is rotating CW
% E0 is the inertial frame located at the origin of the crack center of
% rotation


% (xc,yc) is the position vector of the roller center w.r.t the origin of
% inertial frame (absolute position)
hCam_xc = Rp(1)*cos(theta_vec + beta(1)) + sx_vec.*cos(theta_vec);
hCam_yc = Rp(1)*sin(theta_vec + beta(1)) + sx_vec.*sin(theta_vec);

vCam_xc = Rp(2)*cos(theta_vec + beta(2)) + sy_vec.*cos(theta_vec);
vCam_yc = Rp(2)*sin(theta_vec + beta(2)) + sy_vec.*sin(theta_vec);



nexttile(1,[2,1]);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
ax = gca;
ax.Color = bgColor;
hold on;

axis equal;
% pause(1);

% i=1:10:numel(time_vec);
% rollers_center_plt_array = gobjects(1, numel(i));
% counter = 1;
% for i=1:10:numel(time_vec)
%     rollers_center_plt_array(counter) = plot(xc(i) + r0*cos(theta_vec), yc(i) + r0*sin(theta_vec), ...
%         'HandleVisibility', 'on', 'Color', linkColor);
%     counter = counter + 1;
%     pause(.01);
% end

dsx_dtheta = 1/omega * sx_dot_vec;
d2sx_dtheta2 = (1/omega)^2 * sx_ddot_vec;

dsy_dtheta = 1/omega * sy_dot_vec;
d2sy_dtheta2 = (1/omega)^2 * sy_ddot_vec;

% Note that we parameterize the cam curve by theta.

% derivate with respect to theta
hCam_xc_p = -Rp(1)*sin(theta_vec + beta(1)) + dsx_dtheta.*cos(theta_vec) - sx_vec.*sin(theta_vec);
hCam_yc_p =  Rp(1)*cos(theta_vec + beta(1)) + dsx_dtheta.*sin(theta_vec) + sx_vec.*cos(theta_vec);

% derivate with respect to theta
vCam_xc_p = -Rp(2)*sin(theta_vec + beta(2)) + dsy_dtheta.*cos(theta_vec) - sy_vec.*sin(theta_vec);
vCam_yc_p =  Rp(2)*cos(theta_vec + beta(2)) + dsy_dtheta.*sin(theta_vec) + sy_vec.*cos(theta_vec);



% second derivate with respect to theta
hCam_xc_pp = -Rp(1)*cos(theta_vec + beta(1)) + d2sx_dtheta2.*cos(theta_vec) - ...
    2*dsx_dtheta.*sin(theta_vec) - sx_vec.*cos(theta_vec);
hCam_yc_pp = -Rp(1)*sin(theta_vec + beta(1)) + d2sx_dtheta2.*sin(theta_vec) + ...
    2*dsx_dtheta.*cos(theta_vec) - sx_vec.*sin(theta_vec);

vCam_xc_pp = -Rp(2)*cos(theta_vec + beta(2)) + d2sy_dtheta2.*cos(theta_vec) - ...
    2*dsy_dtheta.*sin(theta_vec) - sy_vec.*cos(theta_vec);
vCam_yc_pp = -Rp(2)*sin(theta_vec + beta(2)) + d2sy_dtheta2.*sin(theta_vec) + ...
    2*dsy_dtheta.*cos(theta_vec) - sy_vec.*sin(theta_vec);



hCam_Delta = 1./((hCam_xc_p.^2 + hCam_yc_p.^2).^(1/2));
hCam_Delta_p = -(hCam_xc_p.^2 + hCam_yc_p.^2).^(-3/2) .* (hCam_xc_p .* hCam_xc_pp + hCam_yc_p .* hCam_yc_pp);


vCam_Delta = 1./((vCam_xc_p.^2 + vCam_yc_p.^2).^(1/2));
vCam_Delta_p = -(vCam_xc_p.^2 + vCam_yc_p.^2).^(-3/2) .* (vCam_xc_p .* vCam_xc_pp + vCam_yc_p .* vCam_yc_pp);


% xCam = xc + rf * yc_p .* Delta;
% yCam = yc - rf * xc_p .* Delta;
%
% xCamp = xc_p + rf * (yc_pp .* Delta + yc_p .* Delta_p);
% yCamp = yc_p - rf * (xc_pp .* Delta + xc_p .* Delta_p);

% Cam Profile
hCam_xCam = hCam_xc - rf(1) * hCam_yc_p .* hCam_Delta;
hCam_yCam = hCam_yc + rf(1) * hCam_xc_p .* hCam_Delta;

vCam_xCam = vCam_xc - rf(2) * vCam_yc_p .* vCam_Delta;
vCam_yCam = vCam_yc + rf(2) * vCam_xc_p .* vCam_Delta;

% rotated by 90 degrees
% vCam_xCam = -(vCam_yc + rf(1) * vCam_xc_p .* vCam_Delta);
% vCam_yCam = vCam_xc - rf(1) * vCam_yc_p .* vCam_Delta;


% rotate_90 = true;
% if rotate_90
%     xx = hCam_xCam;
%     hCam_xCam = -hCam_yCam;
%     hCam_yCam = xx;
% end

hCam_xCamp = hCam_xc_p - rf(1) * (hCam_yc_pp .* hCam_Delta + hCam_yc_p .* hCam_Delta_p);
hCam_yCamp = hCam_yc_p + rf(1) * (hCam_xc_pp .* hCam_Delta + hCam_xc_p .* hCam_Delta_p);

vCam_xCamp = vCam_xc_p - rf(2) * (vCam_yc_pp .* vCam_Delta + vCam_yc_p .* vCam_Delta_p);
vCam_yCamp = vCam_yc_p + rf(2) * (vCam_xc_pp .* vCam_Delta + vCam_xc_p .* vCam_Delta_p);



% now, that the cam x,y coordinates and its derivates are know, we can
% compute instantenous speed of the cam's curve
hCam_S = cumtrapz(theta_vec, sqrt(hCam_xCamp.^2 + hCam_yCamp.^2));
vCam_S = cumtrapz(theta_vec, sqrt(vCam_xCamp.^2 + vCam_yCamp.^2));
% figure;
% subplot(2,1,1);
% plot(time_vec, instantenous_speed_vec);
% subplot(2,1,2);
% plot(time_vec, S);

x0 = -e(2) + r_BA(1) + d(1) + l(1);
y0 = e(1) + r_BA(2) + d(2);

x_tip = x0 + sx_vec;
y_tip = y0 + sy_vec;


%%
% cam1Plt = plot(xCam1, yCam1, 'b-', 'DisplayName', 'Cam1 Profile', 'LineWidth', 1.1);
% cam2Plt = plot(xCam2, yCam2, 'r-', 'DisplayName', 'Cam2 Profile', 'LineWidth', 1.1);
% figure(fg_animation);
% pause(.5);
hCam_camPlt = patch('XData', hCam_xCam, 'YData', hCam_yCam);
vCam_camPlt = patch('XData', vCam_xCam, 'YData', vCam_yCam);
set(hCam_camPlt, 'FaceColor', cam2Color, 'EdgeColor', cam2Color, 'LineStyle', 'none');
set(vCam_camPlt, 'FaceColor', vCamColor, 'EdgeColor', vCamColor, 'LineStyle', 'none');


%
% hold on;
% quiver(xCam2, yCam2, xCam2p, yCam2p, "AutoScaleFactor", .5);


%

% draw base circle with line
hCam_base_circle_plt(1) = plot(rb(1)*cos(theta_vec),rb(1)*sin(theta_vec), 'Color', linkColor, 'DisplayName', 'Base Circle');
hCam_base_circle_plt(2) = plot([0, hCam_xCam(1)], [0, hCam_yCam(1)], 'Color', linkColor);

vCam_base_circle_plt(1) = plot(rb(2)*cos(theta_vec),rb(2)*sin(theta_vec), 'Color', linkColor, 'DisplayName', 'Base Circle');
vCam_base_circle_plt(2) = plot([0, vCam_xCam(1)], [0, vCam_yCam(1)], 'Color', linkColor);


% draw links

hLink = plot([0, l(1)], [0, 0], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', 2, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

% vLink = plot([0, 0], [0, l(2)], 'LineWidth', 2, 'Color', 'k', ...
%     'MarkerIndices', 2, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

vLink(1) = plot([0 0 r_BA(1)  r_BA(1)], [0 l(2), l(2) r_BA(2)], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', 4, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

p = [l(1)/1.3, r_BA(2) + e(1)];
h = .2;
w = .4;

p11 = p - [0,h/2] - [w/2,0];
p12 = p11 + [w,0];
p21 = p11 + [0,h];
p22 = p21 + [w,0];

vLink(2) = plot([0 l(1)/2 p(1)], [l(2), l(2) p11(2)], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4);


vLink(3) = plot([p11(1), p12(1)], [p11(2), p12(2)], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

vLink(4) = plot([p21(1), p22(1)], [p21(2), p22(2)], 'LineWidth', 2, 'Color', 'k', ...
    'MarkerIndices', [], 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 4);




% draw roller
hCam_roller_plt(1) = patch('XData', rf(1)*cos(theta_vec), 'YData', rf(1)*sin(theta_vec), ...
    'FaceColor', linkColor, 'FaceAlpha', 1);

vCam_roller_plt(1) = patch('XData', rf(2)*cos(theta_vec), 'YData', rf(2)*sin(theta_vec), ...
    'FaceColor', linkColor, 'FaceAlpha', 1);

numLines = 6;
for i=1:numLines
    hCam_roller_plt(1+i) = plot([0, rf(1)*cos(i * 2*pi/numLines)], [0 rf(1)*sin(i * 2*pi/numLines)], 'Color', 'k');
    vCam_roller_plt(1+i) = plot([0, rf(2)*cos(i * 2*pi/numLines)], [0 rf(2)*sin(i * 2*pi/numLines)], 'Color', 'k');
end

% tip plot
tip_plt = plot(0,0, 'Color', '#FFD35A', 'LineWidth', 2);
set(tip_plt, 'XData', [], 'YData', []);

% roller_plt(2) = plot([0, -r0], [0 0], 'Color', 'k');


hLinkHTrans = hgtransform;
hCam_HTrans = hgtransform;
hCam_rollerHTrans = hgtransform;

vLinkHTrans = hgtransform;
vCam_HTrans = hgtransform;
vCam_rollerHTrans = hgtransform;

set(hLink, 'Parent', hLinkHTrans, 'Color', linkColor, 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor);
set(vLink, 'Parent', vLinkHTrans, 'Color', linkColor, 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor);

set([hCam_camPlt, hCam_base_circle_plt] , 'Parent', hCam_HTrans);
set([vCam_camPlt, vCam_base_circle_plt] , 'Parent', vCam_HTrans);


set(hCam_roller_plt, 'Parent', hCam_rollerHTrans);
set(vCam_roller_plt, 'Parent', vCam_rollerHTrans);


% draw Cam1 Center
hCam_center = plot(0, 0, 'o', 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor, 'MarkerSize', 8);
vCam_center = plot(0, 0, 'o', 'MarkerFaceColor', linkColor, 'MarkerEdgeColor', linkColor, 'MarkerSize', 8);
%-----------------------------
hLinkHTrans.Matrix = makehgtform('translate', [-e(2)+r_BA(1)+d(1)+sx_vec(1), e(1)+r_BA(2)+d(2)+sy_vec(1), 0]);
vLinkHTrans.Matrix = makehgtform('translate', [-e(2), d(2)+sy_vec(1), 0]);

hCam_HTrans.Matrix = makehgtform('translate', [-e(2)+r_BA(1), d(2)+sy_vec(1)+r_BA(2), 0])*makehgtform('zrotate', -theta_vec(1));
vCam_HTrans.Matrix = makehgtform('zrotate', pi/2 + -theta_vec(1));


hCam_rollerHTrans.Matrix = makehgtform('translate', [-e(2)+r_BA(1)+d(1)+sx_vec(1), e(1)+r_BA(2)+d(2)+sy_vec(1), 0]) * makehgtform('zrotate', -hCam_S(1)/rf(1));
vCam_rollerHTrans.Matrix = makehgtform('translate', [-e(2), d(2)+sy_vec(1), 0]) * makehgtform('zrotate', -vCam_S(1)/rf(2));

axis(ax, [-3.4972    5.4015   -1.6480    4.3844]);
ax.DataAspectRatio = [1 1 1];

xline(0, 'w-');
yline(0, 'w-');

%%
pause(1);
freq = 1/median(diff(time_vec));
t=0; %time initialization
n = 1;
tic;    %start time measuring
% I added the "ishandle" so the program will end in case u closed the figure
while (n <= numel(time_vec))

    hLinkHTrans.Matrix = makehgtform('translate', [-e(2)+r_BA(1)+d(1)+sx_vec(n), e(1)+r_BA(2)+d(2)+sy_vec(n), 0]);
    vLinkHTrans.Matrix = makehgtform('translate', [-e(2), d(2)+sy_vec(n), 0]);

    hCam_HTrans.Matrix = makehgtform('translate', [-e(2)+r_BA(1), d(2)+sy_vec(n)+r_BA(2), 0])*makehgtform('zrotate', -theta_vec(n));
    vCam_HTrans.Matrix = makehgtform('zrotate', pi/2 + -theta_vec(n));

    
    hCam_rollerHTrans.Matrix = makehgtform('translate', [-e(2)+r_BA(1)+d(1)+sx_vec(n), e(1)+r_BA(2)+d(2)+sy_vec(n), 0]) * makehgtform('zrotate', hCam_S(n)/rf(1));
    vCam_rollerHTrans.Matrix = makehgtform('translate', [-e(2), d(2)+sy_vec(n), 0]) * makehgtform('zrotate', vCam_S(n)/rf(2));

    if time_vec(n) <= 5*hh
        set(tip_plt, 'XData', [tip_plt.XData, x_tip(n)], 'YData', [tip_plt.YData, y_tip(n)]);
    end

    set(vline, 'Value', time_vec(n));
    drawnow;
    t = toc; %measuring current time
    n = round(t*freq)+1;
end