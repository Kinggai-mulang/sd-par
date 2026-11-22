function plot_operation_count_bars(dataFile, metric, outputFile)
%PLOT_OPERATION_COUNT_BARS Plot one Experiment 5 event-count panel.

arguments
    dataFile (1, 1) string
    metric (1, 1) string {mustBeMember(metric, ["igt", "swaps"])}
    outputFile (1, 1) string
end

data = readtable(dataFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
if metric == "igt"
    values = [data.CDReductionIGTMean, data.SDTotalIGTMean];
    legendLabels = {'CD reduction', 'SD inheritance + reduction'};
    yLabel = 'Mean integer Gauss transformations';
else
    values = [data.CDReductionSwapMean, data.SDTotalSwapMean];
    legendLabels = {'CD reduction', 'SD inheritance + reduction'};
    yLabel = 'Mean swaps';
end
if any(~isfinite(values), 'all')
    error('plot_operation_count_bars:MissingMetric', ...
        'The selected operation count contains missing or non-finite values.');
end

style = sd_plot_style("two-panel");
fig = figure('Visible', 'off');
ax = axes(fig);
bars = bar(ax, 1:height(data), values, 'grouped', 'BarWidth', 0.78);
bars(1).FaceColor = style.colors.cd;
bars(2).FaceColor = style.colors.sd;
bars(1).EdgeColor = 'none';
bars(2).EdgeColor = 'none';

tickPositions = find(mod(data.Dimension, 10) == 0);
xticks(ax, tickPositions);
xticklabels(ax, compose('%d', data.Dimension(tickPositions)));
xlabel(ax, 'Matrix dimension');
ylabel(ax, yLabel);
legend(ax, legendLabels, 'Location', 'northoutside', ...
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
