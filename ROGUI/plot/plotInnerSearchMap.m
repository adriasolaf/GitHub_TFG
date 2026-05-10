function h = plotInnerSearchMap(ax, innerMap, selectedPoint, metricName)
%plotInnerSearchMap Plot branch-resolved inner DSM samples as points.
%   Multiple values at the same nu are expected when several revolution or
%   Lambert branches are shown. nu is displayed in degrees.
%
%   The color encodes nu_DSM itself as a debugging aid. The vertical axis is
%   the selected metric, usually total delta-v. Samples are not connected by
%   a line because different branches can coexist at the same nu.

    if nargin < 1 || ~isgraphics(ax, 'axes')
        ax = gca;
    end
    if nargin < 2
        innerMap = table();
    end
    if nargin < 3
        selectedPoint = [];
    end
    if nargin < 4 || isempty(metricName)
        metricName = 'totalDv';
    end

    cla(ax);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    h = struct('points', gobjects(0), 'selected', gobjects(0), 'colorbar', gobjects(0));

    if isempty(innerMap) || ~istable(innerMap) || ~any(strcmp('nu', innerMap.Properties.VariableNames))
        showEmptyAxes(ax, 'Inner search map unavailable');
        hold(ax, 'off');
        return;
    end
    if ~any(strcmp(metricName, innerMap.Properties.VariableNames))
        if any(strcmp('dV_DSM', innerMap.Properties.VariableNames))
            metricName = 'dV_DSM';
        else
            metricName = 'totalDv';
        end
    end

    x = rad2deg(innerMap.nu);
    y = innerMap.(metricName);
    if any(strcmp('isFeasible', innerMap.Properties.VariableNames))
        y(~innerMap.isFeasible) = NaN;
    end
    valid = isfinite(x) & isfinite(y);
    if ~any(valid)
        showEmptyAxes(ax, 'Inner search map has no finite values');
        hold(ax, 'off');
        return;
    end

    % Draw individual branch samples only. A line would falsely imply a
    % single-valued continuous curve through different branch families.
    h.points = scatter(ax, x(valid), y(valid), 24, x(valid), 'filled');
    colormap(ax, 'turbo');
    h.colorbar = colorbar(ax);
    ylabel(h.colorbar, '\nu_{DSM} [deg]');
    xlabel(ax, '\nu_{DSM} [deg]');
    ylabel(ax, metricName, 'Interpreter', 'none');
    title(ax, 'Inner branch samples');

    marker = selectedMarker(selectedPoint, innerMap, metricName);
    if ~isempty(marker)
        h.selected(1) = plot(ax, marker(1), marker(2), 'ko', 'LineWidth', 1.8, 'MarkerSize', 12);
        h.selected(2) = plot(ax, marker(1), marker(2), 'rx', 'LineWidth', 2.2, 'MarkerSize', 12);
    end

    hold(ax, 'off');
end

function marker = selectedMarker(selectedPoint, innerMap, metricName)
    % When no explicit point is selected, mark the current best feasible
    % sample. Otherwise mark the point selected by click or auto-selection.
    marker = [];
    if isempty(selectedPoint)
        feasible = isfinite(innerMap.(metricName));
        if any(strcmp('isFeasible', innerMap.Properties.VariableNames))
            feasible = feasible & innerMap.isFeasible;
        end
        if ~any(feasible)
            return;
        end
        rows = innerMap(feasible, :);
        [~, idx] = min(rows.(metricName));
        selectedPoint = rows(idx, :);
    end

    if istable(selectedPoint)
        marker = [rad2deg(selectedPoint.nu(1)), selectedPoint.(metricName)(1)];
    elseif isstruct(selectedPoint) && isfield(selectedPoint, 'nu')
        marker = [rad2deg(selectedPoint.nu), selectedPoint.(metricName)];
    end
end

function showEmptyAxes(ax, message)
    text(ax, 0.5, 0.5, message, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.35 0.35 0.35]);
    xlabel(ax, '\nu_{DSM} [deg]');
    ylabel(ax, '\DeltaV [km/s]');
    title(ax, 'Inner branch samples');
end
