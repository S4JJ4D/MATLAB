%% Example 8.2 of Bertsekas
%
% Romeo and Juliet start dating, but Juliet will be late on any
% date by a random amount X, uniformly distributed over the interval [0, OJ. The
% parameter 0 is unknown and is modeled as the value of a random variable e,
% uniformly distributed between zero and one hour. Assuming that Juliet was late by
% an amount x on their first date, how should Romeo use this information to update
% the distribution of 8?
%
%%
close all;
clear;

% uniform prior
theta_prior_dist = makedist('Uniform', 'lower', 0, 'upper', 1);

% Select a value for theta. remember that theta is constant but unknown.
theta = .8;
% make a uniform dist. for the amount of time Juilet is late for the date
delay_dist = makedist('Uniform', 'lower', 0, 'upper', theta);

x = 0:1e-3:1;
y = pdf(theta_prior_dist, x);
p_plot = plot(x, y, 'b-');
hold on;
xline(theta, 'k--');
grid on;

tobj = title(['Converging To The Truth', newline, 'Dates Counts = 0']);


% hAxes2 = axes(gcf);
% hold on;
% L1 = plot(nan, nan, 'b-');
% L2 = plot(nan, nan, 'r--');
% lgd = legend(hAxes2, [L1, L2], {'Posterior PDF', ['True Value of p = ', num2str(p)]}, ...
%     'Location', 'northeast', 'Box', 'on', 'AutoUpdate', 'off');
% 
% set(hAxes2, 'Visible', 'off');
% set(gcf,'defaultLegendAutoUpdate','off');


% Number of tosses
N = 1000;
x_vec = zeros(N, 1);
x_max = 0;

x_lim_min = .78;
x_lim_max = .82;

xlim([x_lim_min, x_lim_max]);
for i=1:N

    % Arrange a date
    delay = random(delay_dist);
    x_vec(i) = delay;
    x_max = max(x_vec);

    c = 1/integral(@(t)t.^(-i), x_max, 1);
%     theta_vec = 0:1e-4:1;
    theta_vec = x_lim_min:1e-4:x_lim_max;

    pdf_x = @(x) c*(theta_vec.^(-i)).*(x>=x_max) + 0.*(x<x_max);
    expectation_integrand = @(x) c*x.*(x.^(-i));
    EX = integral(expectation_integrand, x_max, 1);

    set(p_plot, 'XData', theta_vec, 'YData', pdf_x(theta_vec));

    tobj.String = ['Converging To The Truth', newline, 'Dates Counts = ', num2str(i), ...
        newline, 'EX = ', num2str(EX)];

    pause(.01);
end
