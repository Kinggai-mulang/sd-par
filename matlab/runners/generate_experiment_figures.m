function files = generate_experiment_figures(mode, experiments, ambiguityOrdering)
%GENERATE_EXPERIMENT_FIGURES Export panels for the selected experiments.

arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    experiments (1, 7) logical = [true true true true false false false]
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

paths = setup_sd_experiment_paths();
dataDirectory = fullfile(paths.numericalOutput, mode);
figureDirectory = fullfile(paths.figureOutput, mode);
if ambiguityOrdering == "random-order"
    dataDirectory = fullfile(dataDirectory, 'random_order');
    figureDirectory = fullfile(figureDirectory, 'random_order');
end
if ~exist(figureDirectory, 'dir'), mkdir(figureDirectory); end

files = strings(0, 1);
if experiments(1)
    dataFile = fullfile(dataDirectory, 'exp01_lambda4_stage_dimension.tsv');
    outputFile = fullfile(figureDirectory, 'exp01_lambda4_stage_dimension.png');
    plot_lambda4_stage_bars(dataFile, "dimension", outputFile);
    files(end+1, 1) = outputFile;
end
if experiments(2)
    dataFile = fullfile(dataDirectory, 'exp02_cd_sd_dimension.tsv');
    outputFile = fullfile(figureDirectory, 'exp02_time_dimension.png');
    plot_cd_sd_bars(dataFile, "dimension", "time", outputFile);
    files(end+1, 1) = outputFile;
    outputFile = fullfile(figureDirectory, 'exp02_flops_dimension.png');
    plot_cd_sd_bars(dataFile, "dimension", "flops", outputFile);
    files(end+1, 1) = outputFile;
end
if experiments(3)
    dataFile = fullfile(dataDirectory, 'exp03_cd_sd_condition.tsv');
    outputFile = fullfile(figureDirectory, 'exp03_time_condition.png');
    plot_cd_sd_bars(dataFile, "condition", "time", outputFile);
    files(end+1, 1) = outputFile;
    outputFile = fullfile(figureDirectory, 'exp03_flops_condition.png');
    plot_cd_sd_bars(dataFile, "condition", "flops", outputFile);
    files(end+1, 1) = outputFile;
end
if experiments(4)
    dataFile = fullfile(dataDirectory, 'exp04_gb.tsv');
    outputFile = fullfile(figureDirectory, 'exp04_gb_time_reduction.png');
    plot_cd_sd_bars(dataFile, "dimension", "time", outputFile);
    files(end+1, 1) = outputFile;
    outputFile = fullfile(figureDirectory, 'exp04_gb_time_search.png');
    plot_cd_sd_bars(dataFile, "dimension", "search", outputFile);
    files(end+1, 1) = outputFile;
end
if experiments(5)
    dataFile = fullfile(dataDirectory, ...
        'exp05_reduction_operations_dimension.tsv');
    outputFile = fullfile(figureDirectory, 'exp05_igt_dimension.png');
    plot_operation_count_bars(dataFile, "igt", outputFile);
    files(end+1, 1) = outputFile;
    outputFile = fullfile(figureDirectory, 'exp05_swaps_dimension.png');
    plot_operation_count_bars(dataFile, "swaps", outputFile);
    files(end+1, 1) = outputFile;
end
if experiments(6)
    dataFile = fullfile(dataDirectory, ...
        'exp06_random_dimension_full_par.tsv');
    outputFile = fullfile(figureDirectory, ...
        'exp06_random_dimension_full_par_ratio.png');
    plot_full_par_time_ratio(dataFile, outputFile);
    files(end+1, 1) = outputFile;
end

if experiments(7)
    dataFile = fullfile(paths.numericalOutput, mode, 'removal_comparison', ...
        'removal_policy_dimension.tsv');
    outputFile = fullfile(paths.figureOutput, mode, 'removal_comparison', ...
        'removal_policy_dimension.png');
    plot_removal_policy_speedup(dataFile, outputFile);
    files(end+1, 1) = outputFile;
end

function plot_full_par_time_ratio(dataFile, outputFile)
data = readtable(dataFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
cdTime = data.CDMeanMs;
sdTime = data.SDMeanMs;
allTimes = [cdTime(:); sdTime(:)];
if any(~isfinite(allTimes)) || any(allTimes <= 0)
    error('generate_experiment_figures:InvalidPARTotalTimes', ...
        'Experiment 6 total PAR times must be finite and positive.');
end
timeRatio = cdTime ./ sdTime;

style = sd_plot_style("single");
fig = figure('Visible', 'off');
ax = axes(fig);
dimensions = data.Dimension;
ratioLine = plot(ax, dimensions, timeRatio, '-o', ...
    'Color', style.colors.speedup, ...
    'LineWidth', style.lineWidth, ...
    'MarkerSize', style.markerSize, ...
    'MarkerFaceColor', style.colors.speedup);
hold(ax, 'on');
yline(ax, 1, '--', 'Color', style.colors.neutral, ...
    'LineWidth', style.referenceLineWidth);
hold(ax, 'off');
tickDimensions = dimensions(mod(dimensions, 10) == 0);
if isempty(tickDimensions)
    tickDimensions = dimensions;
end
xticks(ax, tickDimensions);
xlabel(ax, 'Matrix dimension');
ylabel(ax, 'CD / SD chain runtime');
title(ax, 'Complete deletion-chain runtime ratio');
yValues = [timeRatio; 1];
yPadding = 0.08 * max(range(yValues), 0.1);
ylim(ax, [max(0, min(yValues) - yPadding), max(yValues) + yPadding]);
legend(ax, ratioLine, {'CD / SD'}, ...
    'Location', 'northwest', 'Box', 'off');
apply_sd_plot_style(fig, ax, style);
outputDirectory = fileparts(outputFile);
if strlength(outputDirectory) > 0 && ~exist(outputDirectory, 'dir')
    mkdir(outputDirectory);
end
exportgraphics(fig, outputFile, 'Resolution', style.exportDpi);
close(fig);
end

end
