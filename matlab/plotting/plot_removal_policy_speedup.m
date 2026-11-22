function plot_removal_policy_speedup(dataFile, outputFile)
%PLOT_REMOVAL_POLICY_SPEEDUP Plot three paired single-child removal policies.
arguments
    dataFile (1,1) string
    outputFile (1,1) string
end
data = readtable(dataFile, 'FileType','text','Delimiter','\t', ...
    'VariableNamingRule','preserve');
ratios = [data.RandomSpeedup, data.MaximumVarianceSpeedup, ...
    data.MinimumVarianceSpeedup];
if any(~isfinite(ratios(:))) || any(ratios(:) <= 0)
    error('plot_removal_policy_speedup:InvalidRatios', ...
        'All three removal-policy ratios must be finite and positive.');
end
style = sd_plot_style("single");
fig = figure('Visible','off');
ax = axes(fig);
hold(ax,'on');
% User-requested categorical colors for the three removal policies.
% Reuse the journal palette and keep line style and markers identical.
colors = [style.colors.total; style.colors.speedup; style.colors.cd];
lines = gobjects(3,1);
for policyIndex = 1:3
    lines(policyIndex) = plot(ax, data.Dimension, ratios(:,policyIndex), ...
        '-o', 'Color',colors(policyIndex,:), ...
        'LineWidth',style.lineWidth, 'MarkerSize',style.markerSize, ...
        'MarkerFaceColor',colors(policyIndex,:));
end
yline(ax,1,'--','Color',style.colors.neutral, ...
    'LineWidth',style.referenceLineWidth);
hold(ax,'off');
xticks(ax,data.Dimension(mod(data.Dimension,10)==0));
xPadding = max(1,0.02*range(data.Dimension));
xlim(ax,[min(data.Dimension)-xPadding,max(data.Dimension)+xPadding]);
xlabel(ax,'Matrix dimension');
ylabel(ax,'CD / SD reduction runtime');
title(ax,'Effect of ambiguity removal');
yValues = [ratios(:);1];
padding = 0.08*max(range(yValues),0.1);
ylim(ax,[max(0,min(yValues)-padding),max(yValues)+padding]);
legend(ax,lines,{'Random removal','Maximum variance','Minimum variance'}, ...
    'Location','northwest','Box','off');
apply_sd_plot_style(fig,ax,style);
outputDirectory = fileparts(outputFile);
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
exportgraphics(fig,outputFile,'Resolution',style.exportDpi);
close(fig);
end
