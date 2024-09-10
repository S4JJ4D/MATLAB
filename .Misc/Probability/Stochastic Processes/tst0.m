%%
close all;
clear;

A0 = 1;
w0 = 1;

N = 2e2;

theta_dist = makedist('Uniform', 'lower', 0, 'upper', pi);


theta_samples = random(theta_dist, 1, N);

time_vec = -5:.01:5;
time_mat = repmat(time_vec, N, 1);
y_mat = A0*cos(w0*time_mat + theta_samples.');

% Theoretical Average:
y_bar = -2*A0/pi * sin(w0*time_vec);

subplot(2,1,1);
plot(time_vec, y_mat);
hold on;
plot(time_vec, 1/N * sum(y_mat), 'k-', 'LineWidth', 1.2);
plot(time_vec, y_bar, 'r-', 'LineWidth', 10, 'Color', [1, 0, 0, .3]);

ax2 = subplot(2,1,2);
hold on;
plot(time_vec, y_bar, 'r-', 'LineWidth', 10, 'Color', [1, 0, 0, .3], 'Tag', 'TheoreticalAvgPlot');
avg_plt = plot(nan, nan, 'k-', 'LineWidth', 1.5, 'Tag', 'AvgPlot');
for i=1:N
    plot(time_vec, y_mat(i,:), 'Color', 'g');
    ax2.Children = circshift(ax2.Children, -1);
    set(avg_plt, 'XData', time_vec, 'YData', 1/i * sum(y_mat(1:i, :), 1));

    pause(1e-4);
end




