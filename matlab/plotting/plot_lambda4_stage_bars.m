function plot_lambda4_stage_bars(dataFile, xMode, outputFile)
%PLOT_LAMBDA4_STAGE_BARS Plot LAMBDA 4.0 reduction and search times.

arguments
    dataFile (1, 1) string
    xMode (1, 1) string {mustBeMember(xMode, ["dimension", "condition"])}
    outputFile (1, 1) string
end

data = readtable(dataFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
style = sd_plot_style("two-panel");

fig = figure('Visible', 'off');
ax = axes(fig);
values = [data.ReductionMeanMs, data.SearchMeanMs];
bars = bar(ax, 1:height(data), values, 'grouped', 'BarWidth', 0.72);
bars(1).FaceColor = style.colors.reduction;
bars(2).FaceColor = style.colors.search;
bars(1).EdgeColor = 'none';
bars(2).EdgeColor = 'none';

if xMode == "dimension"
    tickPositions = find(mod(data.Dimension, 5) == 0);
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
ylabel(ax, 'Mean computation time (ms)');
legend(ax, {'Reduction', 'Search'}, 'Location', 'northoutside', ...
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
