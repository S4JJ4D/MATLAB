close all;
clear;

p = .5;
alpha = .9;
% number of samples: time index set size
n = 100;
time_vec = 1:n;

% the mean function
mu = p * (1-alpha.^time_vec)./(1 - alpha);

w = double(rand(1, n) >= p);

subplot(2,1,1);
stem(w);


x = filter(1.0, [1.0 -alpha], w);

subplot(2,1,2);
stem(x, 'DisplayName', 'filtered noise');
hold on;
plot(mu, 'DisplayName', 'mean function');

legend;
