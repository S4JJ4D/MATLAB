%% Story 8.3.3 of Blitzstein
%
% We have a coin that lands Heads with probability p, but we don t know what p is.
% Our goal is to infer the value of p after observing the outcomes of n tosses of 
% the coin. The larger that n is, the more accurately we should be able to estimate p.
%
% There are several ways to go about doing this. One major approach is Bayesian
% inference, which treats all unknown quantities as random variables. In the Bayesian
% approach, we would treat the unknown probability p as a random variable and give
% p a distribution. This is called a prior distribution, and it reflects our uncertainty
% about the true value of p before observing the coin tosses. After the experiment is
% performed and the data are gathered, the prior distribution is updated using Bayes’
% rule; this yields the posterior distribution, which reflects our new
% beliefs about p.
%%
close all;
clear;

prior_dist = makedist('Beta', 'a', 1 ,'b', 1);
posterior_dist = makedist('Beta', 'a', 1 ,'b', 1);

% Select a value for p. remember that p is constant but unknown.
p = .234;
% make a Bernoulli distribution for coin tossing
toss_dist = makedist('Binomial', 'N', 1, 'p', p);

x = 0:1e-3:1;
y = pdf(prior_dist, x);
p_plot = plot(x, y, 'b-');
hold on;
xline(p, 'k--');
grid on;

tobj = title(['Converging To The Truth', newline, 'Toss Counts = 0']);


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
k = 0;
for i=1:N
    % Flip a coin
    toss_result = random(toss_dist);
    k = k + toss_result;
    posterior_dist.a = 1 + k;
    posterior_dist.b = 1 + (i - k);

    y = pdf(posterior_dist, x);
    p_plot.YData = y;

    tobj.String = ['Converging To The Truth', newline, 'Toss Counts = ', num2str(i)];

    pause(.01);
end


