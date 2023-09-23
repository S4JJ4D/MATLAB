function [rest, xp_seq]=tsp(F, t, x0, h, options)
%TSP   p-th Order Taylor Series Numerical Method For Solving ODEs.
%
%   TSP(F,T,X0,H) computes the approximation points {x_n} for the given
%   gradient function F in the time-span specified by t, starting from the
%   initial condition x0. Approxmation values are evaluated at time
%   instants separated by time-step h.
%
%   TSP(..., 'ExactSolution', X) uses the exact solution 'X' to
%   compute global error at each step. Also, the exact solution is plotted
%   if plotting is enabled.
%
%   TSP(..., 'PlotResult', true) plots the resulting approximation
%   values obtained from numerically solving the IVP.
%
%   TSP(..., 'PauseDuration', P) Specifies a the duration of pause (in
%   seconds) within the animation loop of plotting.
%
%   TSP(..., 'PlotInterpolatingCurve', true) Also plots interpolating
%   curves (a polynomial between two successive approximation points
%   directly obtained from TS(P) equations).
%
%   F must be a function handle with a signature of the form @(x,t).
%   In general, F is a matrix-valued function whose columns are the
%   consecutive derivatives of x up to a certain order:
%   F = [f1, f2, f3, ..., fp] where fi is a column vector representing
%   the i'th derivate of the state vector x. The number of columns, p,
%   determines the order of TS method to be used in the function.
%   Therefore, size(F,1) denotes the number of states in the state vector
%   and size(F,2) denotes the order p of the numerical method TS(p).
%   For instance, F = [f1] uses forward Euler method, and F = [f1, f2, f3]
%   uses TS(3) method to compute the approximation to the solution of the
%   IVP.
%
%   
%   REST=TSP(F,T,X0,H) returns the approximation points generated
%   by the numerical method in a tabulated form.
%
%   [REST,XP_SEQ]=TSP(F,T,X0,H) also returns the approximations for the
%   derivative of the IVP solution. XP_SEQ is a 3D matrix where its 1st
%   dimension (along the horizontal axis, to the right) is the time axis, 
%   its 2nd dimension (along the vertical axis towards the bottom) holds
%   the information for the state variables and the 3rd dimension (to the
%   plane) holds the derivatives of higher orders.
%
%   ----------------------------------------------------------------------
%   Example [1]
%
%   F = @(x,t) 2*x*(1-x);
%   x = @(t) 1./(1+4*exp(2*(10-t)));
%   [rest, xp_seq] = tsp(F, [10 13], 1/5, .2, ...
%                        'ExactSolution', x, ...
%                        'PlotResult', true, ...
%                        'PauseDuration', .05, ...
%                        'PlotInterpolatingCurve', true)
%
%   ----------------------------------------------------------------------
%   Example [2]
%
%   F = @(x,t) [[x(2);t-x(1)], ...
%               [t-x(1);1-x(2)], ...
%               [1-x(2);x(1)-t]];
%   x = @(t) [t + cos(t) + sin(t);cos(t) - sin(t) + 1];
%   [rest, ~] = tsp(F, [0 3], [1;2], .1, ...
%                   'ExactSolution', x, ...
%                   'PauseDuration', .05, ...
%                   'PlotResult', true, ...
%                   'PlotInterpolatingCurve', true)
%
%   ----------------------------------------------------------------------
%   Example [3]
%
%   F = @(x,t) [[x(2);x(3);-2*x(1)-x(2)-2*x(3)+(1-2*t)], ...
%               [x(3);-2*x(1)-x(2)-2*x(3)+(1-2*t);4*x(1)+3*x(3)+4*t-4]];
%   x = @(t) [2*sin(t) - t + 1; 2*cos(t) - 1; -2*sin(t)];
%   [rest, ~] = tsp(F, [0 5], [1;1;0], .2, ...
%                   'ExactSolution', x, ...
%                   'PauseDuration', .05, ...
%                   'PlotResult', true, ...
%                   'PlotInterpolatingCurve', true)
%
%   ----------------------------------------------------------------------
%   Example [4]
%
%   F = @(x,t) ...
%       [x(2); -x(1)./(sqrt(x(1).^2 + x(3).^2)).^3;
%        x(4); -x(3)./(sqrt(x(1).^2 + x(3).^2))^3];
%   [rest, ~] = tsp(F, [0 10], [2;0;0;.5], .2, ...
%                   'PauseDuration', .02, ...
%                   'PlotResult', true, ...
%                   'PlotInterpolatingCurve', true)
%   ----------------------------------------------------------------------
%

arguments
    F   (1,1) function_handle {mustBeAFunctionOfNArguments(F, 2, '@(x,t)'), mustBeOfPrescribedForm(F, '@(x,t)')} % a matrix containing gradient functions [f1, f2, ..., fp] where fi is the i'th derivative of x
    t   (2,1) double {mustBeReal, mustBeNonempty, mustHaveNonNegativeLength(t)} % timespan, specified as [t0, tf]
    x0  (:,1) double {mustBeReal, mustBeNonempty, mustBeOfCompatibleSizeWithFcn(x0, F, {'x0', 'F'})} % vector of initial conditions
    h   (1,1) double {mustBeReal, mustBeNonempty, mustBePositive} % time-step

    options.ExactSolution           (:,1) function_handle ...
        {mustBeAFunctionOfNArguments(options.ExactSolution, 1, '@(t)'), ...
        mustBeOfPrescribedForm(options.ExactSolution, '@(t)'), ... % exact analytical function @(t)
        mustBeOfCompatibleSizeWithFcn(x0, options.ExactSolution, {'x0', 'options.ExactSolution'})}
    options.PlotResult              (1,1) logical = false;
    options.PauseDuration           (1,1) {mustBeReal, mustBeNonnegative} = .5;
    options.PlotInterpolatingCurve  (1,1) logical = false;
end

% discretize the time domain:
N = round((t(end)-t(1))/h); % number of time intervals
t0 = t(1);
tf = t0+N*h;
t_seq = t0:h:tf; % note that numel(t_seq) = N+1
% this is because t_seq = {t0, t1, ..., tN}

% define approximation seq and initial condition
n_data_points = N+1;
index_seq = 0:N;
% number of states variables:
[state_count, ts_order] = size(F(x0, 1));

x_seq = zeros(state_count,n_data_points);
x_seq(:,1) = x0;

% x_prime sequence: approximations for the derivatives of the solution curve
xp_seq = zeros(state_count, n_data_points, ts_order);

x_exact_seq = zeros(state_count,n_data_points);
x_exact_seq(:,1) = x0;
GE = zeros(state_count,n_data_points);

is_exact_available = 0;
if isfield(options, 'ExactSolution')
    is_exact_available = 1;
    x = options.ExactSolution;
end

cfs = (((h*ones(1,ts_order)).^(1:ts_order)) ./ factorial(1:ts_order)).';

for i=1:N
    xp_seq(:,i,:) = F(x_seq(:,i), t_seq(i));
    x_seq(:, i+1) = x_seq(:, i) +F(x_seq(:,i), t_seq(i)) * cfs;

    if is_exact_available
        x_exact_seq(:, i+1) = x(t_seq(i+1));
        GE(:, i+1) = x_exact_seq(:, i+1) - x_seq(:, i+1);
    end
end

rest = table(index_seq', t_seq', x_seq', 'VariableNames', {'n', 'tn', 'xn'});
if is_exact_available
    % concat GE if it is available
    rest = [rest, table(GE', 'VariableNames', {'GE'})];
end


%% If PlotResult is enabled
if options.PlotResult
    % open up a figure and start drawing
    % create a new axes
    figure;
    div_list = divisors(state_count);
    if mod(numel(div_list),2) == 0
        idx = [div_list(numel(div_list)/2), div_list(numel(div_list)/2 + 1)];
    else
        idx = [div_list(ceil(numel(div_list)/2)), div_list(ceil(numel(div_list)/2))];
    end

    tlobj = tiledlayout(idx(1), idx(2), 'TileSpacing', 'compact', 'Padding', 'compact');
    if ts_order == 1
        extra_title = ': Forward Euler';
    else
        extra_title = '';
    end
    title(tlobj, sprintf(['TS(%d) Method', extra_title], ts_order));
    subtitle(tlobj, sprintf('Step Size = %.3f', h), 'FontSize', 9);
    for i=1:state_count
        nexttile(i);
        xlabel('Time');
        ylabel('Amplitude');
        title(sprintf('State: x(%d)', i));
        hold on;
        box on;
    end

    if is_exact_available
        hh = .001*h;
        tc = t(1):hh:t(2);
        x_exact = x(tc);
        for i=1:state_count
            ax = nexttile(i);
            plot(ax, tc, x_exact(i,:), 'Color', '#071952', 'Tag', ['Exact Solution: Component ', num2str(i)]);
        end
    end

    if options.PlotInterpolatingCurve
        interpolating_curve_plt = gobjects(1, state_count);
        for i=1:state_count
            ax = nexttile(i);
            interpolating_curve_plt(i) = plot(ax, rest.tn(1), rest.xn(1, i), 'b-', 'Tag', ['Interpolating Curve: Component ', num2str(i)]);
        end
    end


    x_seq_plt = gobjects(1, state_count);
    for i=1:state_count
        ax = nexttile(i);
        x_seq_plt(i) = plot(ax, rest.tn(1), rest.xn(1, i), ...
            'LineStyle', 'none', 'Marker', 'o', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'w', ...
            'MarkerSize', 6, 'Tag', ['Approximation Sequence {x_n}: Component ', num2str(i)]);
    end

    % initial pause
    pause(.5);
    pause_duration = options.PauseDuration;
    N = rest.n(end);
    for i=1:N
        if options.PlotInterpolatingCurve
            hh = .01*h;
            x_int = 0:hh:h;
            y_int = rest.xn(i,:)';
            for j=1:ts_order
                y_int = y_int + 1/factorial(j) * xp_seq(:, i, j) * (x_int.^j);
            end

            for j=1:state_count
                set(interpolating_curve_plt(j), ...
                    'XData', [interpolating_curve_plt(j).XData, rest.tn(i) + x_int], ...
                    'YData', [interpolating_curve_plt(j).YData, y_int(j, :)]);
            end
        end
        for j=1:state_count
            set(x_seq_plt(j), 'XData', [x_seq_plt(j).XData, rest.tn(i+1)], ...
                'YData', [x_seq_plt(j).YData, rest.xn(i+1, j)]);
        end
        if pause_duration ~=0
            pause(pause_duration);
        end
    end
end

end


% ---------------------------------------------------------------------------------
%% Custom validation function

function mustBeAFunctionOfNArguments(f, argcount, fcn_hint)
if nargin(f) ~= argcount
    erridType = 'mustBeAFunctionOfNArguments:NoArgcountTemplateMatch';
    msgType = ['Function handle must accept <strong>', num2str(argcount), '</strong> ' ...
        'input argument(s): ', fcn_hint, '.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustBeOfPrescribedForm(f, form)
f_str = func2str(f);
if ~strcmp(f_str(1:numel(form)), form)
    form_bold = ['<strong>', form, '</strong>'];
    erridType = 'mustBeOfPrescribedForm:NoArgFormTemplateMatch';
    msgType = ['Function handle input argument must be of the form ', ...
        form_bold, '.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustBeOfCompatibleSizeWithFcn(x, f, inputnames)
arguments
    x
    f
    inputnames = {'',''}
end
% checks whether the number of states (row-size of f) matches the number of
% rows in the vector/matrix x.
% remove argument section
f_str = strrep(func2str(f), '@(x,t)', '');
[states_count, ~] = size(str2sym(f_str));
if states_count ~= size(x,1)
    erridType = 'mustBeOfCompatibleSizeWithFcn:InconsistentDimensions';
    msgType = [...
        'Vector ', inputnames{1}, ...
        ' and function handle ', inputnames{2}, ...
        ' are dimensionally inconsistent.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustHaveNonNegativeLength(t)
if t(end) < t(1)
    erridType = 'mustHavePositiveLength:PositiveLength';
    msgType = 'Input vector must be of non-negative length.';
    throwAsCaller(MException(erridType,msgType))
end
end




