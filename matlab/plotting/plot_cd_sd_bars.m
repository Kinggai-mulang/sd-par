function plot_cd_sd_bars(dataFile, xMode, metric, outputFile)
%PLOT_CD_SD_BARS Plot one CD-versus-SD reduction comparison panel.

arguments
    dataFile (1, 1) string
    xMode (1, 1) string {mustBeMember(xMode, ["dimension", "condition"])}
    metric (1, 1) string {mustBeMember(metric, ["time", "flops", "search"])}
    outputFile (1, 1) string
end

data = readtable(dataFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
if metric == "time"
    values = [data.CDMeanMs, data.SDMeanMs];
elseif metric == "flops"
    values = [data.CDMeanFlops, data.SDMeanFlops] / 1e6;
else
    values = [data.CDSearchMeanMs, data.SDSearchMeanMs];
end
if any(~isfinite(values), 'all')
    error('plot_cd_sd_bars:MissingMetric', ...
        'The selected metric contains missing or non-finite values.');
end

style = sd_plot_style("two-panel");
fig = figure('Visible', 'off');
ax = axes(fig);
bars = bar(ax, 1:height(data), values, 'grouped', 'BarWidth', 0.78);
bars(1).FaceColor = style.colors.cd;
bars(2).FaceColor = style.colors.sd;
bars(1).EdgeColor = 'none';
bars(2).EdgeColor = 'none';

if xMode == "dimension"
    dimensionSteps = diff(data.Dimension);
    if height(data) > 10 && all(dimensionSteps == 2)
        % Dense Experiment 2 grid: retain every bar but label decades only.
        tickPositions = find(mod(data.Dimension, 10) == 0);
    else
        % Sparse or irregular grids, including the ten GF dimensions, must
        % display their actual sampled dimensions rather than decade ticks.
        tickPositions = 1:height(data);
    end
    tickLabels = compose('%d', data.Dimension(tickPositions));
    xlabel(ax, 'Matrix dimension');
else
    tickPositions = 1:height(data);
    exponents = round(log2(data.ConditionNumber));
    tickLabels = compose('2^{%d}', exponents);
    xlabel(ax, 'Condition number');
end
xticks(ax, tickPositions);
xticklabels(ax, tickLabels);
if metric == "time" || metric == "search"
    ylabel(ax, 'Mean computation time (ms)');
else
    ylabel(ax, 'Mean FLOPs (\times10^6)');
end
legend(ax, {'CD', 'SD'}, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off');
apply_sd_plot_style(fig, ax, style);
ax.YLim = [0, 1.08 * max(values, [], 'all')];
ax.LooseInset = max(ax.TightInset, [0.02 0.02 0.02 0.02]);

outputDirectory = fileparts(outputFile);
if strlength(outputDirectory) > 0 && ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end
exportgraphics(fig, outputFile, 'Resolution', style.exportDpi);
close(fig);

end
