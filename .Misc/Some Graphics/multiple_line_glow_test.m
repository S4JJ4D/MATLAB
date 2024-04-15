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

xline(0, 'k--');
yline(0, 'k--');

% Create main line

lines = [[2;1], [5;8], [11;7], [12; 0]];
lines_num = size(lines, 2) - 1;

data_struct = repmat(struct(), 1, lines_num);

% p_whole_line = plot(lines(1,:), lines(2,:), 'Color', [1, 0, 0, 0.2], ...
%     'LineWidth', .5, 'Marker', '.', 'MarkerSize', 10);

bars_num = 10; % a measure of smoothness of transition b/w line widths
max_line_width = 4;
min_line_width = .5;
% each bar consists of four points
width_list = linspace(min_line_width, max_line_width, bars_num);

for i=1:lines_num

    d = lines(:,i+1) - lines(:,i);
    l = .065 * norm(d);  % length of each bar
    e = d/norm(d);

    % plot static line
    points = lines(:,i) + d * linspace(0,1,10);
    plot(ax, points(1,:), points(2,:), 'Color', [0, 0, 0, 0.75], ...
    'LineWidth', .5, 'Marker', '.', 'MarkerSize', 5);

    bars_list = gobjects(1, bars_num);
    for j=1:bars_num
        bars_list(j) = plot([0, 0], [0, 0], ...
            'Color', 'r', 'LineWidth', width_list(j));
    end
    supp_bar = plot([0, 0], [0, 0], ...
        'Color', 'r');
    
    data_struct(i).ps = lines(:,i);
    data_struct(i).pe = lines(:,i+1);
    data_struct(i).d = d;
    data_struct(i).l = l;
    data_struct(i).e = e;
    data_struct(i).bars_list = bars_list;
    data_struct(i).supp_bar = supp_bar;

end

%%


axis([-0.7971   14.4480   -2.6092   11.5633]);
% axis equal;

deltas = diff(axis(ax));
delta_x = deltas(1);
delta_y = deltas(3);
delta = sqrt(delta_x^2 + delta_y^2);
scale_val = delta * (1/25);

corners_num = lines_num -1;
for i=1:corners_num
    dif_val = diff(lines(:,i:(i+2)), 1, 2);
    cen_val = lines(:,i+1);
    WayPointCircle(ax, ...
        dif_val(:,1)', dif_val(:,2)', cen_val',...
        'scale', scale_val, 'arrow_count', 3, ...
        'offset', .3, 'arrow_color', '#4b184f');
end


% animate

f = @(t,T) t - floor(t/T)*T;

% frequency
time = 0:.03:10;
for t=time

    lambda0 = f(t, 1);

    for j=1:lines_num

        d = data_struct(j).d;
        l = data_struct(j).l;
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

