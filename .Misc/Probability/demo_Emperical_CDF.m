close all;
fg = figure("Name", 'Demo: Emperical CDF');
fg.Color = [1, 1, 1];

weibull_dist = makedist('Weibull','a',11,'b',7);

xx = 0:.01:15;
yy = cdf(weibull_dist, xx);

sampleSize = 2e3;
sampleSet = random(weibull_dist, 1, sampleSize);

xVals = 0:.05:20;
empericalCDF = zeros(1, numel(xVals));


weibullCDF = cdf(weibull_dist, xVals);
trueCDFPlot = plot(xVals, weibullCDF);
hold on;
grid on;
empericalCDFPlot = plot(0, 0, 'r-', 'LineWidth', 1.5);
set(empericalCDFPlot, 'XData', [], 'YData', []);

rmseText = text(0.7090, 0.7790, {'N = '; 'rmse = '}, 'FontName', 'Source Code Pro', ...
    'EdgeColor', 'k');

ax = gca;
set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontName', 'Source Code Pro', ...
    'FontSize', 9);

title('Estimating CDF by WLLN', 'FontName', 'Source Code Pro Medium');
xlabel('x');
ylabel('F_X(x)');

hAxes2 = axes(gcf);
hold on;
L1 = plot(nan, nan, 'Color', 'b');
L2 = plot(nan, nan, 'Color', 'r', 'LineWidth', 1.5);
lgd = legend(hAxes2, [L1, L2], {'True CDF', 'Emperical CDF'}, ...
    'Location', 'northwest', 'Box', 'on', 'AutoUpdate', 'off');

set(hAxes2, 'Visible', 'off');
set(gcf,'defaultLegendAutoUpdate','off');

% L1 = plot(nan, nan, 'Color', 'b');
% L2 = plot(nan(1,1), 'Color', 'r', 'LineWidth', 1.5);
% lgd = legend([L1, L2], {'True CDF', 'Emperical CDF'}, ...
%     'Location', 'northwest', 'Box', 'on', 'AutoUpdate', 'off');


% lgd = legend({'True CDF', 'Emperical CDF'}, 'Location', 'northwest', 'Box', 'on', ...
%     'AutoUpdate', 'off');

% delete(lgd.AxesListenerList);
% delete(lgd.SelfListenerList);

% setappdata(gca,'LegendColorbarManualSpace',1);
% setappdata(gca,'LegendColorbarReclaimSpace',1);


% Some pre-processing
textStr_N =  "N = " + (1:sampleSize);

pause(.75);
for j=1:sampleSize
    
    currentSampleSet = sampleSet(1:j);
    currentSampleSize = j;

    i = 1;
    for x = xVals
        empericalCDF(i) = 1/currentSampleSize * sum(currentSampleSet < x);
        i = i + 1;
    end
    rmseVal = rmse(empericalCDF, weibullCDF);
    set(rmseText, 'String', {sprintf('N = %d', j); sprintf('rmse = %.5f', rmseVal)});
%     set(rmseText, 'String', [textStr_N(j); "rmse = "+rmseVal]);

    set(empericalCDFPlot, 'XData', xVals, 'YData', empericalCDF);
    pause(1e-3);

end
