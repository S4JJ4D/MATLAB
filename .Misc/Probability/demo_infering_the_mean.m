%% Example 8.3 of Bertsekas

% Inferring the mean using measurements.

% We observe a collection X = (X1, ... Xn) of random variables, with an
% unknown common mean whose value we wish to infer. We assume that given the
% value of the common mean, the Xi are normal and independent, with known variances.
% In a Bayesian approach to this problem, we model the common
% mean as a random variable THETA, with a given prior. For concreteness, we assume a
% normal prior, with known mean x0 and known variance sigma0^2.


%%
close all;
clear;

% normal prior
mu0 = 0;
sigma0 = 1;
theta_prior_dist = makedist('Normal', 'mu', mu0, 'sigma', sigma0);
theta_posterior_dist = makedist('Normal', 'mu', mu0, 'sigma', sigma0);

% Select a value for mu. remember that mu is constant but unknown.
mu_true = 4;

% make a normal dist for X;
% signal is constant and a known value
sigma = 2;
X_dist = makedist('Normal', 'mu', mu_true, 'sigma', sigma);

x_range = -3*sigma0:1e-3:3*sigma0 + mu0;
y = pdf(theta_prior_dist, x_range);
p_plot = plot(x_range, y, 'b-');
hold on;
xline(mu_true, 'k--');
grid on;

tobj = title(['Converging To The Truth', newline, 'Measurements Counts = 0']);


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
N = 1500;
x_vec = zeros(N, 1);
x_max = 0;

x_lim_min = 3.5;
x_lim_max = 4.5;

xlim([x_lim_min, x_lim_max]);
m = mu0;
for i=1:N

    % Make a measurement
    x = random(X_dist);

    m = (i+1-1)/(i+1) * (m + 1/(i+1-1) * x);
    v = (1/(i+1)) * sigma;

    theta_posterior_dist.mu = m;
    theta_posterior_dist.sigma = v;

    x_range = (-5*v:1e-4:5*v) + m;
    y = pdf(theta_posterior_dist, x_range);
    set(p_plot, 'XData', x_range, 'YData', y);


    tobj.String = ['Converging To The Truth', newline, 'Dates Counts = ', num2str(i)];

    pause(.01);
end




