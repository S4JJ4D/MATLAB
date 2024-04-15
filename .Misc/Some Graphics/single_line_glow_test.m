% RUN ME!

if exist('fg', 'var')
    if ishandle(fg)
        delete(fg);
    end
end

fg = figure(...
    'Units', 'normalized', ...
    'Position', [0.6453 0.3417 0.3464 0.5444]);
ax = axes;
hold(ax, 'on');
grid on;
box on;

xline(0, 'r--');
yline(0, 'r--');

% Create main line

ps = [2, 1]';
pe = [5, 8]';

d = pe - ps;
e = d/norm(d);

points_num = 17;
coeffs = linspace(0, 1, points_num);

points = ps + d .* coeffs;
p = plot(points(1,:), points(2,:), 'Color', [1, 0, 0, 0.2], ...
    'LineWidth', .5, 'Marker', '.', 'MarkerSize', 1);

line_ang = atan2(d(2), d(1));

%%
bars_num = 10; % a measure of smoothness of transition b/w line widths
max_line_width = 4;
min_line_width = .5;

l = .065 * norm(d);  % length of each bar

% each bar consists of four points
width_list = linspace(min_line_width, max_line_width, bars_num);

bars_list = gobjects(1, bars_num);
for i=1:bars_num
    bars_list(i) = plot([0, 0], [0, 0], ...
        'Color', 'r', 'LineWidth', width_list(i));
end

supp_bar = plot([0, 0], [0, 0], ...
    'Color', 'r');

starting_point_lambda = .85;
% p0 : starting point on the line
p0 = ps + starting_point_lambda*d;

for i=1:bars_num

    p1 = p0 + l*e;
    if norm(p1 - ps)/norm(d) > 1
        % out of line; compute mod
        alpha = l - norm(pe - p0);
        set(bars_list(i), ...
            'XData', [p0(1), pe(1)], ...
            'YData', [p0(2), pe(2)]);

        set(supp_bar, ...
            'XData', [ps(1), ( ps(1) + alpha*e(1) )], ...
            'YData', [ps(2), ( ps(2) + alpha*e(2) )], ...
            'LineWidth', bars_list(i).LineWidth);


        % update starting point for the next bar
        p0 = ps + alpha * e;
    else
        set(bars_list(i), ...
            'XData', [p0(1), p1(1)], ...
            'YData', [p0(2), p1(2)]);

        % update starting point for the next bar
        p0 = p1;
    end

end



axis([-1.6667    8.7459   -0.0024    9.6776]);

% animate


lambda0 = .2;
% p0 : starting point on the line
p0 = ps + lambda0*d;

f = @(t,T) t - floor(t/T)*T;

% frequency
time = 0:.03:10;
for t=time

    lambda0 = f(t, 1);
    p0 = ps + lambda0*d;

    for i=1:bars_num

        p1 = p0 + l*e;
        if norm(p1 - ps)/norm(d) > 1
            % out of line; compute mod
            alpha = l - norm(pe - p0);
            set(bars_list(i), ...
                'XData', [p0(1), pe(1)], ...
                'YData', [p0(2), pe(2)]);

            set(supp_bar, ...
                'XData', [ps(1), ( ps(1) + alpha*e(1) )], ...
                'YData', [ps(2), ( ps(2) + alpha*e(2) )], ...
                'LineWidth', bars_list(i).LineWidth);


            % update starting point for the next bar
            p0 = ps + alpha * e;
        else
            set(bars_list(i), ...
                'XData', [p0(1), p1(1)], ...
                'YData', [p0(2), p1(2)]);

            % update starting point for the next bar
            p0 = p1;
        end

    end

    drawnow;
%     pause(.05);

end

