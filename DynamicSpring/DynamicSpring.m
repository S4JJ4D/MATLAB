classdef DynamicSpring < handle
    %SPRING represents a spring with varying height
    %   Here, we describe a spring geometric structure
    % by 6 variables (not all are independent) that are related to
    % eachother by the following equations:
    %
    % spring constraint equations
    % 1. s = theta_tot * sqrt(r ^2 + p ^ 2);
    % 2. h = theta_tot * p;
    % 3. theta_tot = 2*pi * k;

    % Changing the radius and length of the spring is allowed. These are
    % considered spring "independent parameters" to be set or changed.
    properties
        name   (1,:) char % spring name
        radius (1,1) double {mustBeReal, mustBeNonnegative} % radius the spring
        pitch  (1,1) double {mustBeReal, mustBeNonnegative} % coil pitch
        turns  (1,1) double {mustBeReal, mustBeNonnegative} % number of coil turns
    end

    % height is considerend the primary variable of the spring to be set by
    % class user. When the height is set externally, spring variables are
    % changed to accomadate the new height. Changes in variables are done
    % in a manner to resemble the strech of compression of the spring when
    % it is mechanically loaded.
    properties (Dependent)
        theta   % total angle through which the spring coil is wound
        length  % total length of the coil
        height  % height of the spring (linear length of the spring)

        hg_tag;
        plt_tag;
    end

    properties(SetAccess = private)
        ax
    end

    methods

        function obj = DynamicSpring(name, radius, pitch, turns)
            %SPRING Construct an instance of this class
            arguments
                name    (1,:) char
                radius
                pitch
                turns
            end
            obj.name = name;
            obj.radius = radius;
            obj.pitch = pitch;
            obj.turns = turns;
        end

        %%
        function SetHeight(obj, height)
            arguments
                obj
                height (1,1) double {mustBeNonnegative, mustBeReal}
            end
            obj.pitch = height/obj.theta;
        end


        %% Plotting
        function PlotSpring(obj, anchorPoint, angle, height, options)
            arguments
                obj
                anchorPoint
                angle
                height
                options.axes (1,1) {mustBeA(options.axes,["matlab.graphics.axis.Axes"," 'matlab.ui.control.UIAxes'"]), mustBeNonDeletedGraphicalObject}
            end

            obj.SetHeight(height);

            if (isempty(obj.ax))
                % no axes are set for plotting
                if (isempty(fieldnames(options)))
                    % no axes is provided by the user either!
                    error('Provide and axes object in the optional arguments to draw the spring');
                end
                % we have a fresh axes to draw objects in
                obj.ax = options.axes;

                hg = hgtransform(obj.ax, 'Tag', obj.hg_tag);
                t = 0:0.01:obj.theta;
                spring_plt = plot(obj.ax, obj.radius * sin(t), obj.pitch*t, 'Tag', obj.plt_tag, 'LineWidth', 1.5);
                set(spring_plt, 'Parent', hg);

                hg.Matrix = makehgtform('translate', [anchorPoint, 0]) * ...
                    makehgtform('zrotate', angle-pi/2);
            else
                % we have axes in our object, but if the user has
                % explicitly supplied an axes object, use that one
                if (~isempty(fieldnames(options)))
                    obj.ax = options.axes;
                else
                    % if the user has not explicitly specified an axes
                    % and we do have an axes in our object, use it:
                end

                % reuse existing stuff if available:

                hg = findobj(obj.ax, 'Tag', obj.hg_tag);
                if (isempty(hg))
                    hg = hgtransform(obj.ax, 'Tag', obj.hg_tag);
                end

                t = 0:0.01:obj.theta;
                spring_plt = findobj(hg, 'Tag', obj.plt_tag);
                if (isempty(spring_plt))
                    spring_plt = plot(obj.ax, obj.radius * sin(t), obj.pitch*t, 'Tag', obj.plt_tag, 'LineWidth', 1.5);
                    set(spring_plt, 'Parent', hg);
                else
                    set(spring_plt, 'XData', obj.radius * sin(t), 'YData', obj.pitch * t);
                end

                hg.Matrix = makehgtform('translate', [anchorPoint, 0]) * ...
                    makehgtform('zrotate', angle-pi/2);
            end

        end


    end

    %% Dependent Methods Getters
    methods
        function theta_val = get.theta(obj)
            theta_val = obj.turns * 2*pi;
        end

        function length_val = get.length(obj)
            length_val = obj.theta * sqrt(obj.radius.^2 + obj.pitch.^2);
        end

        function height_val = get.height(obj)
            height_val = obj.theta * obj.pitch;
        end

        function hg_tag_val = get.hg_tag(obj)
            hg_tag_val = [obj.name, '_hg'];
        end

        function plt_tag_val = get.plt_tag(obj)
            plt_tag_val = [obj.name, '_plt'];
        end

    end
end

%% Custom validation function
function mustBeNonDeletedGraphicalObject(a)
if ~isgraphics(a)
    eidType = 'mustBeNonDeletedGraphicalObject:deletedObject';
    msgType = 'graphical input must be non-deleted';
    throwAsCaller(MException(eidType,msgType))
end
end

