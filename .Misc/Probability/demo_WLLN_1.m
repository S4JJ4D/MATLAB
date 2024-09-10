%% 1
close all;
clear;

% 1. Define the distribution object:
std_uniform_pd = makedist('Uniform');

% 2. Draw a sample of size 1xsample_size with "random" function:

N = 2000;
epsilon = .05;

realizationCount = 10;

for i=1:realizationCount
    %     subplot(3,1,i);
    plt(i) = plot(0, 0, 'LineStyle', '-', 'Color', 'k');
    hold on;
    yline(mean(std_uniform_pd), 'k-');
    ylim([mean(std_uniform_pd) - epsilon, mean(std_uniform_pd) + epsilon]);
end

set(plt, 'XData', [], 'YData', []);
samples = random(std_uniform_pd, realizationCount, N);

Sn = zeros(realizationCount, 1);
Mn = zeros(realizationCount, 1);

for i=1:N
    Sn = Sn + samples(:, i);
    Mn = Sn/i;
    for j=1:realizationCount
        set(plt(j), 'XData', [plt(j).XData, i], 'YData', [plt(j).YData, Mn(j)]);
    end

    if mod(i,100) == 0
        n = i;
        std_normal_pd = makedist('normal', 'mu', 0, 'sigma', sqrt(1/(12*n)));
        xx = -.1:1e-4:.1;
        yy = pdf(std_normal_pd, xx);

        X = [xx;yy;zeros(1,numel(xx));ones(1,numel(xx))];
        X = HTrans('x', n) * HTrans('y', .5) * HRotd('z', -90) * X;

        patch('XData', X(1,:), 'YData', X(2,:), ...
            'FaceColor', '#DAD9EB', 'EdgeColor', '#4A457F', 'FaceAlpha', .3);
    end

    pause(1e-4);
end
