function plot_glowing_lines_2(ax, lines)
%PLOT_GLOWING_LINES Summary of this function goes here
%   Detailed explanation goes here

% lines = [[2;1], [5;8], [11;7], [12; 0]];
% lines = [[2;1], [5;8]];

% Create main line

lines_num = size(lines, 2) - 1;

data_struct = repmat(struct(), 1, lines_num);

p_whole_line = plot(ax, lines(1,:), lines(2,:), 'Color', [1, 0, 0, 0.2], ...
    'LineWidth', .5, 'Marker', '.', 'MarkerSize', 10);

bars_num = 10; % a measure of smoothness of transition b/w line widths
max_line_width = 3;
min_line_width = .5;
l = .5;  % length of each bar
% each bar consists of four points
width_list = linspace(min_line_width, max_line_width, bars_num);

for i=1:lines_num

    d = lines(:,i+1) - lines(:,i);
    e = d/norm(d);

    % plot static line
    points = lines(:,i) + d * linspace(0,1,10);
    plot(ax, points(1,:), points(2,:), 'Color', [1, 0, 0, 0.2], ...
    'LineWidth', .5, 'Marker', '.', 'MarkerSize', 5);


    bars_list = gobjects(1, bars_num);
    for j=1:bars_num
        bars_list(j) = plot(ax, [0, 0], [0, 0], ...
            'Color', 'r', 'LineWidth', width_list(j));
    end
    supp_bar = plot(ax, [0, 0], [0, 0], ...
        'Color', 'r');
    
    data_struct(i).ps = lines(:,i);
    data_struct(i).pe = lines(:,i+1);
    data_struct(i).d = d;
    data_struct(i).e = e;
    data_struct(i).bars_list = bars_list;
    data_struct(i).supp_bar = supp_bar;

end

%%


axis([-0.7971   14.4480   -2.6092   11.5633]);

% animate

f = @(t,T) t - floor(t/T)*T;

% frequency
time = 0:.03:10;
for t=time

    lambda0 = f(t, 1);

    for j=1:lines_num

        d = data_struct(j).d;
        e = data_struct(j).e;
        ps = data_struct(j).ps;
        pe = data_struct(j).pe;
        bars_list = data_struct(j).bars_list;
        supp_bar = data_struct(j).supp_bar;


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
    end

    drawnow;
    %     pause(.05);

end


end

