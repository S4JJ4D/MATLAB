%%
clear;
close all;

% k: number of coil turns
% b is changed to adapt h
k = 10;
t = 0:.01:(k*2*pi);
r = .1;
b = 1;

% total length
s = k * 2*pi * sqrt(r.^2 + b.^2);

f = @(t) [r*cos(t); r*sin(t); b*t];

data = f(t);

subplot(1,3,1);
sprintPlt = plot3(data(1,:), data(2,:), data(3,:));
axis([-r r -r r 0 1.1]);
ax.DataAspectRatio = [1 1 1];
view(-40, 6.5);

subplot(1,3,2);
plt2 = plot(data(1,:), data(3,:), 'LineWidth', 1.1);
axis([-r r 0 1.1]);
ax.DataAspectRatio = [1 1 1];


subplot(1,3,3);
plt3 = plot(data(2,:), data(3,:), 'LineWidth', 1.1);
axis([-r r 0 1.1]);
ax.DataAspectRatio = [1 1 1];




ht = 0:.1:20;
hv = 5 + 5*sin(ht);

for h=hv
    % as h changes, ...
    lambda = h/s;
    b = r * (lambda)/sqrt(1-lambda^2);

    set(sprintPlt, 'ZData', b*t);
    set([plt2, plt3], 'YData', b*t);
    pause(.02);

end




