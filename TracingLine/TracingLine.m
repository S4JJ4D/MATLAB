classdef TracingLine < handle
    %TracingLine represents a dynamic curve with limited buffer size to
    %hold trajectory points. curve length is thus limited as its motion
    %evolves in time.

    % Example:
    % 1. Create an instance

    properties
        name        (1,:) char
        bufferSize  (1,1) double {mustBeReal, mustBeNonnegative}

        xBuffer     (1,:) double
        yBuffer     (1,:) double

        % An emtpy ax variable is interpreted as an uset ax variable
        ax {mustBeA(ax,["matlab.graphics.axis.Axes", "matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
    end


    properties (Dependent)
        linePltTag; 
    end

    methods
        function obj = TracingLine(name, bufferSize, options)
            %TracingLine Construct an instance of this class
            arguments
                name        (1,:) char
                bufferSize  (1,1) double {mustBeReal, mustBePositive}
                options.xBufferInit (1,:) double = []
                options.yBufferInit (1,:) double = []
                options.axes {mustBeA(options.axes,["matlab.graphics.axis.Axes","matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
            end

            if (any(string(TracingLine.NamesVessel()) == string(name)))
                error(['Duplicate name! Use another name. Names already in use are: ', newline, sprintf('\t->\t%s\n', string(TracingLine.NamesVessel()))])
            else
                obj.name = name;
                TracingLine.NamesVessel('add', string(name));
            end

            obj.bufferSize = bufferSize;

            if(length(options.xBufferInit) > bufferSize)
                warning('Supplied XData array is larger than the buffer size. Truncating leading entries in XData array to fit the buffer.');
                obj.xBuffer = options.xBufferInit(end-bufferSize+1:end);
            elseif(length(options.xBufferInit) < bufferSize)
                warning('Supplied XData array is smaller than the buffer size. Prepending zero entries to XData array to fit the buffer.');
                n = length(options.xBufferInit);
                obj.xBuffer = [zeros(1, bufferSize-n), options.xBufferInit];
            else
                obj.xBuffer = options.xBufferInit;
            end

            if(length(options.yBufferInit) > bufferSize)
                warning('Supplied YData array is larger than the buffer size. Truncating leading entries in YData array to fit the buffer.')
                obj.yBuffer = options.yBufferInit(end-bufferSize+1:end);
            elseif (length(options.yBufferInit) < bufferSize)
                warning('Supplied YData array is smaller than the buffer size. Prepending zero entries to YData array to fit the buffer.');
                n = length(options.yBufferInit);
                obj.yBuffer = [zeros(1, bufferSize-n), options.yBufferInit];
            else
                obj.yBuffer = options.yBufferInit;
            end

            obj.ax = options.axes;
        end
    end

    methods
        % remove the name from the NamesVessel upon deletion
        function delete(obj)
            TracingLine.NamesVessel("delete", string(obj.name));
        end
    end

    % ------------------------------------------------------------------
    methods
        function AddPoint(obj, xData, yData, updatePlot)
            arguments
                obj
                xData (1,:) double = []
                yData (1,:) double = []
                updatePlot (1,1) logical = false
            end
            xN = length(xData);
            obj.xBuffer = circshift(obj.xBuffer, -xN);
            obj.xBuffer(end-xN+1:end) = xData;

            yN = length(yData);
            obj.yBuffer = circshift(obj.yBuffer, -yN);
            obj.yBuffer(end-yN+1:end) = yData;

            if(updatePlot)
                obj.PlotLine();
            end

        end

        function PlotLine(obj, options)
            arguments
                obj
                options.Color     (1,:) char = 'red'
                options.LineWidth (1,1) double = 1.2
                options.HeadMarker (1,1) logical = false
                % Additional arguments could be added for extra plotting options
            end
            if (isempty(obj.ax))
                error('''obj.ax'' member is unset. Set the ax member before plotting');
            end

            % reuse existing graphic objects if available:
            tracingLinePlt = findobj(obj.ax, 'Tag', obj.linePltTag);
            if (isempty(tracingLinePlt))
                % create a fresh plt object if none already exists
                tracingLinePlt = plot(obj.ax, obj.xBuffer, obj.yBuffer, 'Tag', obj.linePltTag, 'LineWidth', options.LineWidth, 'Color', options.Color);
                if (options.HeadMarker)
                    set(tracingLinePlt, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', options.Color, 'MarkerIndices', obj.bufferSize)
                end
            else
                set(tracingLinePlt, 'XData', obj.xBuffer, 'YData', obj.yBuffer);
            end

        end

    end

    % ------------------------------------------------------------------
    % get dependent method
    methods
        function linePltTag_val = get.linePltTag(obj)
            linePltTag_val = [obj.name, '_plt'];
        end
    end

    % ------------------------------------------------------------------
    % tracks names to avoid conflict
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
    % ------------------------------------------------------------------

end
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
        msgType = 'An input of type double for ''axes'' argument must be empty.';
        throwAsCaller(MException(eidType,msgType))
    end
end
end

