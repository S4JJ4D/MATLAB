classdef DynamicSpring < handle
    %DYNAMICSPRING    Graphical spring object with adaptive length.
    %   DYNAMICSPRING(name, 'Radius', r, 'Pitch', p, 'Turns', t) constructs and
    %   initializes a spring object with the given name and parameters.
    %
    %   DYNAMICSPRING(..., 'Configuration', config) positions the spring
    %   object at the given anchor point specified by config(1:2) and
    %   orients it on the plane at the given angle config(3)
    %
    %   DYNAMICSPRING(..., 'Axes', ax) attaches a specific axes object to
    %   the spring so that plots are made within this particular axes
    %   object.
    %
    %   DYNAMICSPRING(..., 'VisualForm', vf) plots the spring object with the
    %   speficed visual form. vf is either "simplified" or "detailed"
    %
    %   Geometric structure of a spring (state) is described by six variables
    %   (not all are independent) that are related via a set of constraint equations
    %   described as:
    %   - s = t * sqrt(r^2 + p^2);
    %   - h = t * p;
    %   - t = 2*pi * k;
    %
    %   Where the state variables are
    %   1. s:      total length of the spring coil
    %   2. t:      total angle through which the spring coil is wound
    %   3. r:      radius of the spring
    %   4. p:      coil pitch
    %   5. h:      linear length/height of the spring
    %   6. k:      number of coild turns
    %
    %   Variables <strong>h</strong>, <strong>p</strong>, and <strong>r</strong>, are chosen as <strong>independent variables</strong>.
    %   These are the variables that can directly be set or modifiedby the user
    %   of the class.
    %
    %   The graphical model of the spring is obtained from the mathematical
    %   descrpition of a rising sinusoidal curve in the 3D space which is
    %   described by the following map: (<strong>t</strong> is the independent
    %   variable)
    %   
    %   t -> (A0*cos(w*t), B0*sin(w*t), C0*t) \in R³
    %
    %   Changing the radius, pitch, and number of turns of a spring is permissible
    %   through the class interface. These are considered spring "defining parameters"
    %   to be set or changed by class user.
    %
    %   Examples:
    %
    %     - Define a spring:
    %       sp1 = DynamicSpring('spring1', 'Radius', .25, 'Pitch', .03, 'Turns', 18, 'Axes', ax);
    %
    %     - Plot the spring when it is anchored at [1, 2], oriented at 30 degrees
    %       w.r.t the horizontal line and is stretched to a given length of 5:
    %       sp1.PlotSpring(5, 'Configuration', [1, 2, pi/6]);
    %
    %     - Plot the spring with "axled" visual form:
    %       sp1.visual_form = "axled";
    %
    %     - Plot the spring with gray color:
    %       sp1.plotting_options.SpringBodyColor = [.6 .6 .6]


    properties
        % An emtpy ax variable is interpreted as an uset obj.ax member
        ax {mustBeA(ax,["matlab.graphics.axis.Axes", "matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
        plotting_options = ...
            struct(...
            'LineWidth', 2, 'SpringBodyColor', 'k', ...
            'RearEyeOuterColor', '#525252', 'RearEyeOuterColorAlpha', 1.0,...
            'FrontEyeOuterColor', '#525252', 'FrontEyeOuterColorAlpha', 1.0,...
            'RearEyeInnerColor', 'none', 'RearEyeInnerColorAlpha', 1.0, ...
            'FrontEyeInnerColor', 'none', 'FrontEyeInnerColorAlpha', 1.0, ...
            'RearClevisColor', 'w', 'FrontClevisColor', 'w');
        visual_form {mustBeMember(visual_form, ["simplified", "detailed", "axled"])} = "detailed"
    end

    properties (Dependent)
        theta              % total angle through which the spring coil is wound
        coil_length        % total length of the coil
        spring_length      % linear length of the spring part of the graphical spring

        is_plotted         % a logical value representing the plot state of the spring
        hg_tag;            % a tag attached to hgtransform object whose child is spring plot
        plt_tag;           % a tag attached to a spring plot object whose parent is a hgtransform
    end

    properties (SetAccess = private)
        % A name must be supplied to distinguish spring plot objects
        name            (1,:) char {mustBeNonempty} = '-'            % the unique name assigned to the object
        radius          (1,1) double {mustBeReal, mustBeNonnegative} % radius of the spring
        pitch           (1,1) double {mustBeReal, mustBeNonnegative} % coil pitch
        turns           (1,1) double {mustBeReal, mustBeNonnegative} % number of coil turns
        total_length    (1,1) double {mustBeReal}                    % total linear length of the graphical spring
        configuration   (1,3) double {mustBeReal}                    % the planar configuration (position and orientation) of the spring
    end

    % ---------------------------------------------------------------------------
    methods
        function obj = DynamicSpring(name, options)
            %DynamicSpring Constructor
            arguments
                name                            (1,:) char
                options.Radius                  (1,1) double {mustBeReal, mustBeNonnegative} = 0.25
                options.Pitch                   (1,1) double {mustBeReal, mustBeNonnegative} = 0.03
                options.Turns                   (1,1) double {mustBeReal, mustBeNonnegative} = 12
                options.Configuration           (1,3) double {mustBeReal}                    = [0 0 0];
                options.Axes {mustBeA(options.Axes,["matlab.graphics.axis.Axes","matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
                options.VisualForm {mustBeMember(options.VisualForm, ["simplified", "detailed", "axled"])} = "detailed"
            end

            % Check for name conflicts:
            if (any(string(DynamicSpring.NamesVessel()) == string(name)))
                error(['Duplicate name! Use another name. Names already in use are: ', newline, sprintf('\t->\t%s\n', string(DynamicSpring.NamesVessel()))])
            else
                DynamicSpring.NamesVessel('add', string(name));
            end

            obj.name = name;
            obj.radius = options.Radius;
            obj.pitch = options.Pitch;
            obj.turns = options.Turns;
            obj.configuration = options.Configuration;
            obj.ax = options.Axes;
            obj.visual_form = options.VisualForm;
        end

        % Remove the name from the NamesVessel upon deletion
        function delete(obj)
            DynamicSpring.NamesVessel("delete", string(obj.name));
        end

        % ---------------------------------------------------------------------------
        %% Plotting
        function PlotSpring(obj, total_length, options)
            arguments
                obj
                total_length            % total length of the graphical spring
                options.Configuration   % [anchor_point, angle]
            end

            if (isempty(obj.ax))
                error('''obj.ax'' member is unset. Set the ''obj.ax'' member before plotting ');
            elseif ~isgraphics(obj.ax)
                obj.ax = [];
                error('''obj.ax'' is a deleted graphics object. Either set the ''obj.ax'' member before plotting ');
            end

            % Adjusting spring shape for the given linear length specified
            % by the user
            obj.total_length = total_length;
            obj.SetSpringLength(total_length);

            % Reuse existing objects if available:
            hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
            if (isempty(hg))
                % create a fresh hg object if none exists
                hg = hgtransform(obj.ax, 'Tag', obj.hg_tag);
            end

            springPlotGp = findobj(hg.Children, 'flat', 'Tag', obj.plt_tag);
            if (isempty(springPlotGp))
                % create a fresh plt object if none exists
                springPlotGp = obj.PlotSpringGraphics();
                springPlotGp.Tag = obj.plt_tag;
                set(springPlotGp, 'Parent', hg);
            else
                obj.UpdateSpringLength(springPlotGp)
            end
            
            % The following could be done using input.parser.
            if any(cellfun(@(x)strcmpi(x,'Configuration'), fieldnames(options)))
                % If a Config is provided by the user
                obj.configuration = options.Configuration;
            end

            % 'makehgtform' is slow. Directly assigning the ht matrix:
            alpha = obj.configuration(3);
            hg.Matrix = ...
                [cos(alpha), -sin(alpha), 0, obj.configuration(1);
                 sin(alpha),  cos(alpha), 0, obj.configuration(2);
                          0,           0, 1,                    0;
                          0,           0, 0,                    1];

        end
    end

    % ---------------------------------------------------------------------------
    %% Customized Setter Methods
    % To prevent calls during the construction phase, a set of customized
    % setter methods is introduced here. Built-in set.method functions are
    % called during the object initialization phase which might be
    % undesirable in some cases.
    methods
        % ----------------------- setName
        function setName(obj, new_name)

            % Check for name conflicts:
            if (any(string(DynamicSpring.NamesVessel()) == string(new_name)))
                error(['Duplicate name! Use another name. Names already in use are: ', newline, sprintf('\t->\t%s\n', string(DynamicSpring.NamesVessel()))])
            else
                % Delete the previous name from the vessel
                DynamicSpring.NamesVessel("delete", string(obj.name));
                % Add the new name into the vessel
                DynamicSpring.NamesVessel('add', string(new_name));
            end

            % graphical plot objects that depend on the name must be
            % updated:
            if (obj.is_plotted)
                % 1. Updating hgtform object
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                hg.Tag = [new_name, '_hg'];
                % 2. Updating plt object
                plt = findobj(hg.Children, 'flat', 'Tag', obj.plt_tag);
                plt.Tag = [new_name, '_plt'];
            end
            obj.name = new_name;
        end


        % ----------------------- setConfiguration
        function setConfiguration(obj, new_config)
            obj.configuration = new_config;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                hg.Matrix = makehgtform('translate', [new_config(1:2), 0]) * ...
                    makehgtform('zrotate', new_config(3));
            end
        end


        % ----------------------- setRadius
        function setRadius(obj, new_radius)
            obj.radius = new_radius;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                % delete existing plot
                delete(findobj(hg.Children, 'flat', 'Tag', obj.plt_tag));
                % create a fresh plt object if none already exists
                springPlotGp = obj.PlotSpringGraphics();
                springPlotGp.Tag = obj.plt_tag;
                set(springPlotGp, 'Parent', hg);
            end
        end

        % ----------------------- setTotalLength
        function setTotalLength(obj, new_total_length)

            % Adjusting spring shape for the given height
            obj.total_length = new_total_length;
            obj.SetSpringLength(new_total_length);

            if obj.is_plotted
                springPlotGp = findobj(obj.ax.Children, 'Tag', obj.plt_tag);
                obj.UpdateSpringLength(springPlotGp);
            end
        end

        % ----------------------- setPitch
        function setPitch(obj, new_pitch)

            obj.pitch = new_pitch;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                % delete existing plot
                delete(findobj(hg.Children, 'flat', 'Tag', obj.plt_tag));
                % create a fresh plt object if none already exists
                springPlotGp = obj.PlotSpringGraphics();
                springPlotGp.Tag = obj.plt_tag;
                set(springPlotGp, 'Parent', hg);
            end
        end


        % ----------------------- setVisualForm
        function setVisualForm(obj, new_visual_form)
            arguments
                obj
                new_visual_form {mustBeMember(new_visual_form, ["simplified", "detailed", "axled"])}
            end
            obj.visual_form = new_visual_form;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                % delete existing plot
                delete(findobj(hg.Children, 'flat', 'Tag', obj.plt_tag));
                % create a fresh plt object if none already exists
                springPlotGp = obj.PlotSpringGraphics();
                springPlotGp.Tag = obj.plt_tag;
                set(springPlotGp, 'Parent', hg);
            end
        end

    end

    % ---------------------------------------------------------------------------
    %% Setter Methods
    methods
        % The primary purpose of this set method
        % is to set 'hold' to 'on' once and for all.
        function set.ax(obj, new_ax)
            arguments
                obj
                new_ax {mustBeA(new_ax,["matlab.graphics.axis.Axes", "matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
            end
            obj.ax = new_ax;
            if (~isempty(new_ax))
                hold(new_ax, 'on');
            end
        end

        function set.plotting_options(obj, new_plotting_options)
            % Warnings can be safely ignored because this method is not
            % called during object constrcution (due to the default value
            % set for the property).
            obj.plotting_options = new_plotting_options;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                % delete existing plot
                delete(findobj(hg.Children, 'flat', 'Tag', obj.plt_tag));
                % create a fresh plt object if none already exists
                springPlotGp = obj.PlotSpringGraphics();
                springPlotGp.Tag = obj.plt_tag;
                set(springPlotGp, 'Parent', hg);
            end
        end
    end

    % ---------------------------------------------------------------------------
    %% Dependent Methods Getters
    methods
        function theta_val = get.theta(obj)
            theta_val = obj.turns * 2*pi;
        end

        function length_val = get.coil_length(obj)
            length_val = obj.theta * sqrt(obj.radius.^2 + obj.pitch.^2);
        end

        function spring_length_val = get.spring_length(obj)
            spring_length_val = obj.theta * obj.pitch;
        end

        function hg_tag_val = get.hg_tag(obj)
            hg_tag_val = [obj.name, '_hg'];
        end

        function plt_tag_val = get.plt_tag(obj)
            plt_tag_val = [obj.name, '_plt'];
        end

        function is_plotted_val = get.is_plotted(obj)
            is_plotted_val = true;
            if isempty(obj.ax)
                is_plotted_val = false;
            elseif ~isgraphics(obj.ax)
                obj.ax = [];
                is_plotted_val = false;
            else
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                if (isempty(hg))
                    is_plotted_val = false;
                end
            end
        end
    end

    % ---------------------------------------------------------------------------
    % tracks names to avoid name conflicts
    methods (Static=true, Access=private)
        function out = NamesVessel(options)
            arguments
                options.add (1,1) string = ""
                options.delete (1,1) string = ""
            end
            persistent NamesList;
            if (strlength(options.add) ~= 0)
                NamesList = [NamesList, options.add];
            end
            if (strlength(options.delete) ~= 0)
                NamesList(string(NamesList) == options.delete) = [];
            end
            if (strlength(options.add) == 0 && strlength(options.delete) == 0)
                out = NamesList;
            end
        end
    end
    % ---------------------------------------------------------------------------
    %%
    methods (Access=private)
        % spring length is considerend the primary variable of the spring to be set by
        % class user. When the height is set externally, spring variables are
        % changed to accomadate the new height. Changes in variables are done
        % in a manner to resemble the strech or compression of the spring when
        % it is mechanically loaded.
        function SetSpringLength(obj, new_spring_length)
            arguments
                obj
                new_spring_length (1,1) double {mustBeNonnegative, mustBeReal}
            end
            obj.pitch = new_spring_length/obj.theta;
        end
    end

    % ---------------------------------------------------------------------------
    %%
    methods (Access=private)

        function springPlotGp = PlotSpringGraphics(obj)

            x_end = obj.total_length;
            size = 20/9 * obj.radius;
            r_eye_out = size/2;
            r_eye_in = r_eye_out/3;
            radius_eye_vec = [r_eye_in, r_eye_out];
            line_width = obj.plotting_options.LineWidth;

            if obj.visual_form == "simplified"
                t = 0:0.01:obj.theta;
                springPlotGp = plot(obj.ax, obj.pitch*t, -obj.radius*sin(t), ...
                    'Tag', 'hg_graphical_spring', 'LineWidth', line_width, 'Color', obj.plotting_options.SpringBodyColor);

            elseif obj.visual_form == "detailed"

                % ------------------------ Group all plots together
                springPlotGp = hggroup();
                springPlotGp.Tag = 'hg_graphical_spring';


                % ------------------------ Case Mounting Joint
                case_patch(1) = PlotAnnulusPatch(obj.ax, [0, 0], radius_eye_vec,...
                    'FaceColor', obj.plotting_options.RearEyeOuterColor,...
                    'FaceAlpha', obj.plotting_options.RearEyeOuterColorAlpha,...
                    'InnerCircleFaceColor', obj.plotting_options.RearEyeInnerColor,...
                    'InnerCircleFaceAlpha', obj.plotting_options.RearEyeInnerColorAlpha,...
                    'N', 60, ...
                    'LineWidth', line_width, 'Tag', 'Case Mounting Joint 1');


                % ------------------------ Case Mounting Clevis 1
                t = linspace(0, pi, 30);
                a = 1.1*r_eye_out;
                b = .3*a;
                x = [-b*sin(t), 1.2*r_eye_out, 1.2*r_eye_out, 0] + r_eye_out;
                y = [a*cos(t), -a, a, a];
                case_patch(2) = patch(obj.ax, 'XData', x, 'YData', y,...
                    'FaceColor', obj.plotting_options.RearClevisColor,...
                    'LineWidth', line_width, ...
                    'Tag', 'Case Mounting Clevis 1');


                % ------------------------ Spring
                t_i = r_eye_out/obj.pitch;
                t_f = (x_end - r_eye_out)/obj.pitch;
                t = t_i:0.05:t_f;
                spring_plt = plot(obj.ax, obj.pitch*t, -obj.radius*sin(t), ...
                    'LineWidth', line_width, 'Color', obj.plotting_options.SpringBodyColor, ...
                    'Tag', 'Spring');
                obj.ax.Children = circshift(obj.ax.Children, -1);


                % ------------------------ Case Mounting Joint 2
                case_patch(3) = PlotAnnulusPatch(obj.ax, [x_end, 0], radius_eye_vec,...
                    'FaceColor', obj.plotting_options.FrontEyeOuterColor,...
                    'FaceAlpha', obj.plotting_options.FrontEyeOuterColorAlpha,...
                    'InnerCircleFaceColor', obj.plotting_options.FrontEyeInnerColor,...
                    'InnerCircleFaceAlpha', obj.plotting_options.FrontEyeInnerColorAlpha,...
                    'N', 60, ...
                    'LineWidth', line_width, 'Tag', 'Case Mounting Joint 2');


                % ------------------------ Case Mounting Clevis 2
                t = linspace(0, pi, 30);
                a = 1.1*r_eye_out;
                b = .3*a;
                x = -[-b*sin(t), 1.2*r_eye_out, 1.2*r_eye_out, 0] + x_end - r_eye_out;
                y = [a*cos(t), -a, a, a];
                case_patch(4) = patch(obj.ax, 'XData', x, 'YData', y,...
                    'FaceColor', obj.plotting_options.FrontClevisColor,...
                    'LineWidth', line_width, ...
                    'Tag', 'Case Mounting Clevis 2');


                set([case_patch, spring_plt], 'Parent', springPlotGp);
                obj.ax.Children(1).Children = circshift(obj.ax.Children(1).Children, -1);

            else

                % ------------------------ Group all plots together
                springPlotGp = hggroup();
                springPlotGp.Tag = 'hg_graphical_spring';


                % ------------------------ Case Mounting Joint
                case_patch(1) = PlotAnnulusPatch(obj.ax, [0, 0], radius_eye_vec,...
                    'FaceColor', obj.plotting_options.RearEyeOuterColor,...
                    'FaceAlpha', obj.plotting_options.RearEyeOuterColorAlpha,...
                    'InnerCircleFaceColor', obj.plotting_options.RearEyeInnerColor,...
                    'InnerCircleFaceAlpha', obj.plotting_options.RearEyeInnerColorAlpha,...
                    'N', 60, ...
                    'LineWidth', line_width, 'Tag', 'Case Mounting Joint 1');


                % ------------------------ Case Block
                x = [0, r_eye_out, r_eye_out, 0] + 1.2*r_eye_out;
                y = [-1.5*r_eye_in, -1.5*r_eye_in, 1.5*r_eye_in, 1.5*r_eye_in];
                case_patch(2) = patch(obj.ax, 'XData', x, 'YData', y, ...
                    'LineJoin', 'round', 'FaceColor', 'w', 'LineWidth', line_width, ...
                    'Tag', 'Case Block');


                % ------------------------ Case Mounting Clevis 1
                t = linspace(0, pi, 30);
                a = 1.1*r_eye_out;
                b = .3*a;
                x = [-b*sin(t), 1.2*r_eye_out, 1.2*r_eye_out, 0] + 2*r_eye_out;
                y = [a*cos(t), -a, a, a];
                case_patch(3) = patch(obj.ax, 'XData', x, 'YData', y,...
                    'FaceColor', obj.plotting_options.RearClevisColor,...
                    'LineWidth', line_width, ...
                    'Tag', 'Case Mounting Clevis 1');


                % ------------------------ Spring
                t_i = 2*r_eye_out/obj.pitch;
                t_f = (x_end - 1.5*r_eye_out)/obj.pitch;
                t = t_i:0.05:t_f;
                spring_plt = plot(obj.ax, obj.pitch*t, -obj.radius*sin(t), ...
                    'LineWidth', line_width, 'Color', obj.plotting_options.SpringBodyColor, ...
                    'Tag', 'Spring');
                obj.ax.Children = circshift(obj.ax.Children, -1);


                % ------------------------ Case Mounting Joint 2
                case_patch(4) = PlotAnnulusPatch(obj.ax, [x_end, 0], radius_eye_vec,...
                    'FaceColor', obj.plotting_options.FrontEyeOuterColor,...
                    'FaceAlpha', obj.plotting_options.FrontEyeOuterColorAlpha,...
                    'InnerCircleFaceColor', obj.plotting_options.FrontEyeInnerColor,...
                    'InnerCircleFaceAlpha', obj.plotting_options.FrontEyeInnerColorAlpha,...
                    'N', 60, ...
                    'LineWidth', line_width, 'Tag', 'Case Mounting Joint 2');


                % ------------------------ Case Mounting Clevis 2
                t = linspace(0, pi, 30);
                a = 1.1*r_eye_out;
                b = .3*a;
                x = -[-b*sin(t), 1.2*r_eye_out, 1.2*r_eye_out, 0] + x_end - 1.5*r_eye_out;
                y = [a*cos(t), -a, a, a];
                case_patch(5) = patch(obj.ax, 'XData', x, 'YData', y,...
                    'FaceColor', obj.plotting_options.FrontClevisColor,...
                    'LineWidth', line_width, ...
                    'Tag', 'Case Mounting Clevis 2');


                % ------------------------ Longitudinal Shaft
                x = [.5*(r_eye_in+r_eye_out), x_end-.5*(r_eye_in+r_eye_out), ...
                    x_end-.5*(r_eye_in+r_eye_out), .5*(r_eye_in+r_eye_out)];
                y = [-r_eye_in, -r_eye_in, r_eye_in, r_eye_in];
                case_patch(6) = patch(obj.ax, 'XData', x, 'YData', y, ...
                    'FaceColor', [0.28,0.22,0.22], 'LineWidth', line_width, ...
                    'Tag', 'Shaft');
                obj.ax.Children = circshift(obj.ax.Children, -1);


                set([spring_plt, case_patch], 'Parent', springPlotGp);
                obj.ax.Children(1).Children = circshift(obj.ax.Children(1).Children, -1);
                obj.ax.Children(1).Children([6,7]) = obj.ax.Children(1).Children([7,6]);

            end

        end

        % ------------------------------------------------------------------------------
        % ------------------------------------------------------------------------------

        function [] = UpdateSpringLength(obj, springPlotGp)


            x_end = obj.total_length;

            if obj.visual_form == "simplified"
                t = 0:0.01:obj.theta;
                set(springPlotGp, 'XData', obj.pitch*t, 'YData', -obj.radius*sin(t));

            elseif obj.visual_form == "detailed"

                % ------------------------ Updating Spring
                spring_plt = findobj(springPlotGp.Children, 'flat', 'Tag', 'Spring');
                t_i = 10/9 * obj.radius/obj.pitch;
                t_f = (x_end - 10/9 * obj.radius)/obj.pitch;
                t = t_i:0.05:t_f;
                spring_plt.XData = obj.pitch*t;
                spring_plt.YData = -obj.radius*sin(t);


                % ------------------------ Updating Case Mounting Joint 2
                case_joint_2 = findobj(springPlotGp.Children, 'flat', 'Tag', 'Case Mounting Joint 2');
                % x_current = mean(case_joint_2.Children(1).XData);
                data = case_joint_2.Children(1).XData;
                x_current = .5 * (data(end/2) + data(1));
                case_joint_2.Children(1).XData = case_joint_2.Children(1).XData + (x_end - x_current);
                case_joint_2.Children(2).XData = case_joint_2.Children(2).XData + (x_end - x_current);


                % ------------------------ Updating Case Mounting Clevis 2
                case_collor_2 = findobj(springPlotGp.Children, 'flat', 'Tag', 'Case Mounting Clevis 2');
                set(case_collor_2, 'XData', case_collor_2.XData + (x_end - x_current));

            else
                
                r_eye_out = 10/9 * obj.radius;

                % ------------------------ Updating Spring
                spring_plt = findobj(springPlotGp.Children, 'flat', 'Tag', 'Spring');
                t_i = 2*r_eye_out/obj.pitch;
                t_f = (x_end - 1.5*r_eye_out)/obj.pitch;
                t = t_i:0.05:t_f;
                spring_plt.XData = obj.pitch*t;
                spring_plt.YData = -obj.radius*sin(t);


                % ------------------------ Updating Case Mounting Joint 2
                case_joint_2 = findobj(springPlotGp.Children, 'flat', 'Tag', 'Case Mounting Joint 2');
                % x_current = mean(case_joint_2.Children(1).XData);
                data = case_joint_2.Children(1).XData;
                x_current = .5 * (data(end/2) + data(1));
                case_joint_2.Children(1).XData = case_joint_2.Children(1).XData + (x_end - x_current);
                case_joint_2.Children(2).XData = case_joint_2.Children(2).XData + (x_end - x_current);


                % ------------------------ Updating Case Mounting Clevis 2
                case_collor_2 = findobj(springPlotGp.Children, 'flat', 'Tag', 'Case Mounting Clevis 2');
                set(case_collor_2, 'XData', case_collor_2.XData + (x_end - x_current));


                % ------------------------ Longitudinal Shaft
                shaft = findobj(springPlotGp.Children, 'flat', 'Tag', 'Shaft');
                xdata = shaft.XData;
                set(shaft, 'XData', ...
                    [xdata(1),...
                     xdata(2) + (x_end - x_current), ...
                     xdata(3) + (x_end - x_current), ...
                     xdata(4)]);

            end
        end

    end

end
% ---------------------------------------------------------------------------------
%% Helper Drawing Functions
function annulusPatch = PlotAnnulusPatch(ax, center, radius, options)
arguments
    ax {mustBeA(ax,["matlab.graphics.axis.Axes","matlab.ui.control.UIAxes"]), mustBeNonDeletedGraphicalObject}
    center (1,2) double {mustBeReal}
    radius (1,2) double {mustBeReal}
    options.N (1,1) {mustBeInteger, mustBeNonnegative} = 30;

    options.LineWidth = 1;
    options.LineStyle = '-';

    options.FaceColor = 'k';
    options.FaceAlpha = 1;

    options.EdgeColor = 'k';
    options.EdgeAlpha = 'flat';

    options.InnerCircleFaceColor = 'w';
    options.InnerCircleFaceAlpha = 0.0;

    options.Tag = '';

end

N = options.N;
t = linspace(0, 2*pi, N);
x = [radius(1)*cos(t) + center(1), radius(2)*cos(t) + center(1)];
y = [radius(1)*sin(t) + center(2), radius(2)*sin(t) + center(2)];


alpha_data = ones(2*N,1);
alpha_data([N, end]) = [0 0]';

outer_circle = patch(ax, ...
    'XData', x, 'YData', y, 'FaceVertexAlphaData', alpha_data, ...
    'LineWidth', options.LineWidth, 'LineStyle', options.LineStyle, ...
    'FaceColor', options.FaceColor, 'FaceAlpha', options.FaceAlpha,...
    'EdgeColor', options.EdgeColor, 'EdgeAlpha', options.EdgeAlpha, ...
    'Tag', [options.Tag, ' outer circle']);


x = radius(1)*cos(t) + center(1);
y = radius(1)*sin(t) + center(2);

inner_circle = patch(ax, ...
    'XData', x, 'YData', y, ...
    'LineStyle', 'none', ...
    'FaceColor', options.InnerCircleFaceColor, ...
    'FaceAlpha', options.InnerCircleFaceAlpha,...
    'Tag', [options.Tag, ' inner circle']);


annulusPatch = hggroup('Tag', options.Tag);
set([outer_circle, inner_circle], 'Parent', annulusPatch);

end

% ---------------------------------------------------------------------------------
%% Custom validation function
function mustBeNonDeletedGraphicalObject(a)
if ~isgraphics(a)
    eidType = 'mustBeNonDeletedGraphicalObject:deletedObject';
    msgType = 'graphical input must be non-deleted.';
    throwAsCaller(MException(eidType,msgType))
end
end

function mustBeEmptyIfDouble(a)
if isa(a, 'double')
    if (~isempty(a))
        eidType = 'mustBeEmptyIfDouble:NonEmptyDouble';
        msgType = ['An input of type <strong>double</strong> for ' ...
            '<strong>axes</strong> argument must be an empty vector.'];
        throwAsCaller(MException(eidType,msgType))
    end
end
end

