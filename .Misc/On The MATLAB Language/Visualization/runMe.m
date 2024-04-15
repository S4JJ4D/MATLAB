

% wingSpan = .175;
% textWingspan = .2;
% offset = .05;
% group1 = hggroup;
% 
% %Visualizing relative forces produced by motors:
% F_A = quiver3(0,-wingSpan,0,0,0,.2,...
%     'LineWidth',2.5,'Color','red','Parent',group1);
% F_B = quiver3(0, wingSpan,0,0,0,.2,...
%     'LineWidth',2.5,'Color','red','Parent',group1);
% F_C = quiver3(wingSpan,0,0,0,0,.2,...
%     'LineWidth',2.5,'Color','red','Parent',group1);
% F_D = quiver3(-wingSpan,0,0,0,0,.2,...
%     'LineWidth',2.5,'Color','red','Parent',group1);

sim('Control_1D');
quadPlot_1D(logsout);