function style = sd_plot_style(layout)
%SD_PLOT_STYLE Journal figure settings for the sequential-decorrelation study.
%   STYLE = SD_PLOT_STYLE(LAYOUT) returns the physical size, typography,
%   line settings, and semantic color palette for one exported image.
%   LAYOUT is 'single', 'two-panel', or 'three-panel' and refers to the
%   image's intended placement in LaTeX.

arguments
    layout (1, 1) string {mustBeMember(layout, ["single", "two-panel", "three-panel"])}
end

style.layout = layout;
style.fontName = 'Times New Roman';
style.exportDpi = 600;
style.backgroundColor = [1 1 1];
style.axesColor = [0.12 0.12 0.12];
style.gridColor = [0.72 0.72 0.72];
style.gridAlpha = 0.45;
style.axesLineWidth = 0.8;
style.referenceLineWidth = 1.0;

% The palette is deliberately small and keeps the same meaning throughout
% the manuscript. Values follow the established colors in the legacy plots.
style.colors.sd = [0.0000 0.6200 0.7800];
style.colors.cd = [0.4940 0.1840 0.5560];
style.colors.reduction = [0.8500 0.3250 0.0980];
style.colors.search = [0.4660 0.6740 0.1880];
style.colors.total = [0.0000 0.4470 0.7410];
style.colors.speedup = [0.4660 0.6740 0.1880];
style.colors.neutral = [0.25 0.25 0.25];

switch layout
    case "single"
        style.figureWidthIn = 3.45;
        style.figureHeightIn = 2.65;
        style.fontSizeTitle = 9.5;
        style.fontSizeLabel = 9.0;
        style.fontSizeTick = 8.0;
        style.fontSizeLegend = 7.5;
        style.lineWidth = 1.25;
        style.markerSize = 4.5;
    case "two-panel"
        style.figureWidthIn = 3.30;
        style.figureHeightIn = 2.45;
        style.fontSizeTitle = 9.0;
        style.fontSizeLabel = 8.5;
        style.fontSizeTick = 7.5;
        style.fontSizeLegend = 7.0;
        style.lineWidth = 1.15;
        style.markerSize = 4.0;
    case "three-panel"
        style.figureWidthIn = 2.15;
        style.figureHeightIn = 1.85;
        style.fontSizeTitle = 8.0;
        style.fontSizeLabel = 7.5;
        style.fontSizeTick = 6.5;
        style.fontSizeLegend = 6.3;
        style.lineWidth = 1.0;
        style.markerSize = 3.2;
end

end
