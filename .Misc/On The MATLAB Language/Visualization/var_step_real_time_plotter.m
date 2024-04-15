t=0; %time initialization
u = @(p) p>0;
n = 1;  %counter
tic;    %start time measuring
% I added the "ishandle" so the program ends in case the figure is closed
while (n <= numel(time))
    
    p_obj.XData = [p_obj.XData, tau(n)];
    p_obj.YData = [p_obj.YData, omega(n)];
    %     p_obj.XData = tau(n);
    %     p_obj.YData = omega(n);
    drawnow;  %updates the display
    title_obj.String = sprintf('time = %.2f',t);
    t = toc; %measuring current time
    n = find(abs(time - t) == min(abs(time - t)))+1;
    %     n = round(t*freq)+1;
end
title_obj.String = sprintf('time = %.2f',time(end));

end