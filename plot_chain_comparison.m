function plot_chain_comparison(mode)
%PLOT_CHAIN_COMPARISON Plot generated GB and random-LDLT sequence results.
arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
end

directory = fileparts(mfilename('fullpath'));
addpath(fullfile(directory, 'matlab'));
config = sd_experiment_config();
paths = setup_sd_experiment_paths();
dataDirectory = fullfile(paths.numericalOutput, mode);
outputDirectory = fullfile(paths.figureOutput, mode);
if config.defaultOrdering == "random-order"
    dataDirectory = fullfile(dataDirectory, 'random_order');
    outputDirectory = fullfile(outputDirectory, 'random_order');
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
filenames = ["exp06_gb_full_par.tsv", "exp06_random_dimension_full_par.tsv"];
models = ["GB", "Random-LDLT"];
style = sd_plot_style("single");
colors = [style.colors.cd; style.colors.sd];
for metric = ["runtime", "search_share"]
    fig = figure('Visible', 'off');
    ax = axes(fig);
    hold(ax, 'on');
    for model = 1:2
        data = readtable(fullfile(dataDirectory, filenames(model)), ...
            'FileType', 'text', 'Delimiter', '\t');
        values = data.CDMeanMs ./ data.SDMeanMs;
        if metric == "search_share"
            values = 100 * data.CDSearchChainMeanMs ./ data.CDMeanMs;
        end
        plot(ax, data.Dimension, values, '-o', 'Color', colors(model, :), ...
            'LineWidth', style.lineWidth, 'MarkerSize', style.markerSize);
    end
    xlabel(ax, 'Starting ambiguity dimension', 'FontName', style.fontName, ...
        'FontSize', style.fontSizeLabel);
    if metric == "runtime"
        ylabel(ax, 'Total-runtime ratio (CD/SD)', 'FontName', style.fontName, ...
            'FontSize', style.fontSizeLabel);
        yline(ax, 1, '--', 'Color', style.colors.neutral, ...
            'HandleVisibility', 'off');
    else
        ylabel(ax, 'Search share of CD runtime (%)', 'FontName', style.fontName, ...
            'FontSize', style.fontSizeLabel);
    end
    apply_sd_plot_style(fig, ax, style);
    legend(ax, models, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
        'Box', 'off', 'FontName', style.fontName, 'FontSize', style.fontSizeLegend);
    exportgraphics(fig, fullfile(outputDirectory, "chain_" + metric + ".png"), ...
        'Resolution', style.exportDpi);
    close(fig);
end
end
