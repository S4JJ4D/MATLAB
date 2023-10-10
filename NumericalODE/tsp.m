function [rest, xp_seq]=tsp(F, t, x0, h, tol, options)
%TSP   p-th Order Taylor Series Numerical Method For Solving ODEs.
%
%   TSP(F,T,X0,H) computes the approximation points {x_n} for the given
%   gradient function F in the time-span specified by t, starting from the
%   initial condition x0. Approxmation values are evaluated at time
%   instants separated by time-step h.
%
%   TSP(F,T,X0,[],tol) uses adaptive step-size to compute approximate
%   solutions for the given IVP. LTE at each time-point is maintained lower
%   than the tolerance specified via TOL argument. If adaptive step-size is
%   activated, the last column of F is used for LTE estimation.
%
%   TSP(..., 'ExactSolution', X) uses the exact solution 'X' to
%   compute global error at each step. Also, the exact solution is plotted
%   if plotting is enabled.
%
%   TSP(..., 'PlotResult', true) plots the resulting approximation
%   values obtained from numerically solving the IVP.
%
%   TSP(..., 'PauseDuration', P) Specifies the duration of pause (in
%   seconds) within the animation loop of plotting.
%
%   TSP(..., 'PlotInterpolatingCurve', true) Also plots interpolating
%   curves (a polynomial between two successive approximation points
%   directly obtained from TS(P) equations).
%
%   TSP(..., 'ReportLTE', true) Also report local truncation error in the
%   output. This option is only effective when an exact solution is
%   available.
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
%   LTE at each time point t_{n+1} is obtained by subtracting the exact
%   value of the solution at that time z(t_{n+1}) from its approximation
%   which is obtained using a truncated version of Taylor expansion.
%   Realize that LTE does NOT use the approximation values {xn} obtained
%   from the numerical method. LTE only accounts for a portion of GE.
%
%
%                                  ^
%                                  z(t+h)
%                      ╭──────────────────────────────╮
%                      │                              │
%
%        z(t+h)   =   z(t) + h*z'(t) + h^2/2! * z''(t)   +  LTE(t+h)
%
%      │        │      │                              │    │        │
%      ╰────────╯      ╰──────────────────────────────╯    ╰────────╯
%     exact solution        estimate of the exact        estimation error
%     at t=t_{n+1}          solution at t=t_{n+1}        committed for this
%                                                        approximation
%
%
%
%
%                                 estimate╭──◎  △
%                            ╭────────────╯     │
%                            │                  │  LTE(n+1)
%                            │                  │
%                     ╭──────╯       ╭───────◎  ▽
%                    ●┴──────────────╯exact
%
%               ─────┃────────────────────────┃───────────────────────▶
%                    ┃                        ┃
%     Time Grid:    t_{n}                    t_{n+1}
%
%
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

%   TODO: use a preallocation technique in adaptive step-size method

arguments
    F   (1,1) function_handle {mustBeAFunctionOfNArguments(F, 2, '@(x,t)'), mustBeOfPrescribedSignature(F, '@(x,t)')} % a matrix containing gradient functions [f1, f2, ..., fp] where fi is the i'th derivative of x
    t   (2,1) double {mustBeReal, mustBeNonempty, mustHaveNonNegativeLength(t)} % timespan, specified as [t0, tf]
    x0  (:,1) double {mustBeReal, mustBeNonempty, mustBeOfCompatibleSizeWithFcn(x0, F, {'x0', 'F'})} % vector of initial conditions
    h         double {mustBeReal, mustBePositive, mustBeScalarOrEmpty} = [] % time-step
    tol       double {mustBeReal, mustBePositive, mustBeScalarOrEmpty} = [] % tolerance

    options.ExactSolution           (:,1) function_handle ...
        {mustBeAFunctionOfNArguments(options.ExactSolution, 1, '@(t)'), ...
        mustBeOfPrescribedSignature(options.ExactSolution, '@(t)'), ... % exact analytical solution @(t)
        mustBeOfCompatibleSizeWithFcn(x0, options.ExactSolution, {'x0', 'options.ExactSolution'})}
    options.PauseDuration           (1,1) {mustBeReal, mustBeNonempty, mustBeNonnegative} = .5;
    options.PlotResult              (1,1) logical = false;
    options.PlotInterpolatingCurve  (1,1) logical = false;
    options.ReportLTE               (1,1) logical = false;

end

%% Input Validation

is_adaptive_step_size_mode = false;
tabchar = char(9);
if isempty(h) && isempty(tol)
    msg = ...
        ['Either a step-size (h) or a tolerance (tol) must be specified.', ...
        newline, tabchar, ...
        'For a fixed step-size approach, use the syntax:    <strong>TSP(F,T,X0,H)</strong>', ...
        newline, tabchar, ...
        'For a variable step-size approach, use the syntax: <strong>TSP(F,T,X0,[],TOL)</strong>'];
    error(msg);
elseif isempty(h) && ~isempty(tol)
    is_adaptive_step_size_mode = true;
elseif ~isempty(h) && isempty(tol)
else
    msg = ...
        ['Both a step-size (h) and a tolerance (tol) cannot be simultaneously specified.', ...
        newline, tabchar, ...
        'For a fixed step-size approach, use the syntax:    <strong>TSP(F,T,X0,H)</strong>', ...
        newline, tabchar, ...
        'For a variable step-size approach, use the syntax: <strong>TSP(F,T,X0,[],TOL)</strong>'];
    error(msg);
end

% number of states variables and the order, p, of TS(p) method:
[state_count, ts_order] = size(F(x0, 1));

% If adaptive step-size mode is activated, the last column of F is used to
% estimate LTE at each time instant. As a result, the first (n-1) columns
% of F is used to advance the numerical solutions.
if is_adaptive_step_size_mode && ts_order < 2
    msg = 'When adaptive step-size mode is enabled, at least two derivative functions must be supplied: size(F,2) > 1';
    error(msg);
elseif is_adaptive_step_size_mode && ts_order >= 2
    % the last column is reserved for LTE estimation
    ts_order = ts_order - 1;
end


%% Obtaining Approximate Solutions

is_exact_available = 0;
if isfield(options, 'ExactSolution')
    is_exact_available = 1;
    x = options.ExactSolution;
end
report_lte = 0;
if options.ReportLTE && is_exact_available
    report_lte = 1;
elseif options.ReportLTE && ~is_exact_available
    warning('LTE cannot be computed unless an exact solution is availabe.');
end

if ~is_adaptive_step_size_mode

    % discretize the time domain
    N = round((t(end)-t(1))/h); % number of time intervals
    t0 = t(1);
    tf = t0+N*h;
    t_seq = t0:h:tf; % note that numel(t_seq) = N+1
    % this is because t_seq = {t0, t1, ..., tN}

    % define approximation seq and initial condition
    n_data_points = N+1;
    index_seq = 0:N;

    x_seq = zeros(state_count,n_data_points);
    x_seq(:,1) = x0;

    % x_prime sequence: approximations for the derivatives of the solution curve
    xp_seq = zeros(state_count, n_data_points, ts_order);
    % populate the first element:
    xp_seq(:,1,:) = F(x_seq(:,1), t_seq(1));

    x_exact_seq = zeros(state_count,n_data_points);
    x_exact_seq(:,1) = x0;
    GE = zeros(state_count,n_data_points);
    LTE = zeros(state_count, n_data_points);

    % coefficients in the Taylor series expansion upto the order specified by
    % ts_order
    cfs = (((h*ones(1,ts_order)).^(1:ts_order)) ./ factorial(1:ts_order)).';

    % Obtaining Approximate Solutions
    % --------------------------------------------------------
    % Main Loop: States are evolved within the following loop
    for i=1:N
        x_seq(:, i+1) = x_seq(:, i) + F(x_seq(:,i), t_seq(i)) * cfs;
        xp_seq(:,i+1,:) = F(x_seq(:,i+1), t_seq(i+1));

        if is_exact_available
            x_exact_seq(:, i+1) = x(t_seq(i+1));
            GE(:, i+1) = x_exact_seq(:, i+1) - x_seq(:, i+1);
        end

        if report_lte
            LTE(:, i+1) = x_exact_seq(:, i+1) - (x_exact_seq(:, i) +F(x_exact_seq(:,i), t_seq(i)) * cfs);
        end
    end
    % --------------------------------------------------------

    rest = table(index_seq', t_seq', x_seq', 'VariableNames', {'n', 'tn', 'xn'});
    if is_exact_available
        % concat GE if it is available
        rest = [rest, table(GE.', 'VariableNames', {'GE'})];
    end
    if report_lte
        rest = [rest, table(LTE.', 'VariableNames', {'LTE'})];
    end


else

    p = ts_order;
    t0 = t(1);
    tf = t(end);

    t_seq = [];
    x_seq = [];
    xp_seq= [];
    h_seq = [];
    x_exact_seq = [];
    GE = [];
    LTE = [];
    LTE_estimate = [];

    lstcol = [zeros(p,1);1];
    h0 = ...
        min(abs((tol*factorial(p+1)) ./ (F(x0, t0)*lstcol)).^(1/(p+1)));

    t_seq = [t_seq, t0];
    x_seq = [x_seq, x0];
    xp_seq = zeros(state_count, 1, p);

    der = F(x0, t0);
    xp_seq(:,1,:) = reshape(der(:,1:end-1), [state_count,1,p]);
    h_seq = [h_seq, h0];
    x_exact_seq = [x_exact_seq, x0];
    GE = [GE, zeros(state_count,1)];
    LTE = [LTE, zeros(state_count, 1)];
    LTE_estimate = [LTE_estimate, zeros(state_count, 1)];

    step_rejection_seq = 0;
    % Obtaining Approximate Solutions
    % --------------------------------------------------------
    % Main Loop: States are evolved within the following loop
    while t_seq(end) < tf
        
        % coefficients in the Taylor series expansion upto the order specified by
        % ts_order
        cfs = (((h_seq(end)*ones(1,p)).^(1:p)) ./ factorial(1:p)).';
        der = F(x_seq(:, end), t_seq(end));
        x_next = x_seq(:, end)+der(:,1:end-1)*cfs;

        LTE_estimate_next = 1/factorial(p+1) * h_seq(end)^(p+1) * der(:,end);
        h_next = h_seq(end) * (abs(tol/max(abs(LTE_estimate_next)))^(1/(p+1)));
        t_next = t_seq(end) + h_seq(end);

        if max(abs(LTE_estimate_next)) > 1.1*tol
            h_seq(end) = h_next;
            step_rejection_seq(end) = 1;
            continue;
        else
            % The current step is accepted.
            % "Taking" the step ...
            t_seq = [t_seq, t_next];
            h_seq = [h_seq, h_next];
            x_seq = [x_seq, x_next];
            LTE_estimate = [LTE_estimate, LTE_estimate_next];

            step_rejection_seq = [step_rejection_seq, 0];

            der = F(x_seq(:, end), t_seq(end));
            xp_seq = [xp_seq, reshape(der(:,1:end-1), [state_count,1,p])];

            if is_exact_available
                x_exact_seq = [x_exact_seq, x(t_seq(end))];
                GE = [GE, x_exact_seq(:, end) - x_seq(:, end)];
            end

            if report_lte
                cfs = (((h_seq(end-1)*ones(1,p)).^(1:p)) ./ factorial(1:p)).';
                der = F(x_exact_seq(:, end-1), t_seq(end-1));
                LTE = [LTE, x_exact_seq(:, end) - (x_exact_seq(:, end-1) + der(:,1:end-1)*cfs)];
            end


        end
    end
    % --------------------------------------------------------

    % Construct the output table by reporting relevant results
    index_seq = 0:numel(t_seq)-1;
    rest = table(index_seq', h_seq', t_seq', x_seq', LTE_estimate.', 'VariableNames', {'n', 'hn', 'tn', 'xn', 'LTE (Estimate)'});
    
    if report_lte
        % concat LTE if available
        rest = [rest, table(LTE.', 'VariableNames', {'LTE'})];
    end

    if is_exact_available
        % concat GE if available
        rest = [rest, table(GE.', 'VariableNames', {'GE'})];
    end

    % Add comment column
    comment = repmat("", size(rest,1), 1);
    % Indicates that the corresponding step was rejected because computed 
    % LTE exceeded the given tolerance.
    comment(step_rejection_seq ~= 0) = "Rejected";
    rest = [rest, table(comment, 'VariableNames', {'Comment'})];

end

% ---------------------------------------------------------------------------------
%% Plot Results

% if plotting is enabled
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
    if isempty(h)
        subtitle(tlobj, 'Adaptive Step-Size', 'FontSize', 9);
        h=.1; % for subsequent plotting purposes
    else
        subtitle(tlobj, sprintf('Step Size = %.3f', h), 'FontSize', 9);
        N = rest.n(end);
        h_seq = h*ones(1,N); % for subsequent plotting purposes
    end
    
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
            hh = .02*h_seq(i);
            x_int = 0:hh:h_seq(i);
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

function mustBeOfPrescribedSignature(f, form)
f_str = func2str(f);
if ~strcmp(f_str(1:numel(form)), form)
    form_bold = ['<strong>', form, '</strong>'];
    erridType = 'mustBeOfPrescribedSignature:NoArgFormTemplateMatch';
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
