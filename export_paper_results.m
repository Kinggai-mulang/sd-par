%% Export the generated tables and figures to paper_results.
scriptDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDirectory, 'matlab'));
config = sd_experiment_config();
paths = setup_sd_experiment_paths();
sourceData = fullfile(paths.numericalOutput, config.experimentMode);
sourceFigures = fullfile(paths.figureOutput, config.experimentMode);
if config.defaultOrdering == "random-order"
    sourceData = fullfile(sourceData, 'random_order');
    sourceFigures = fullfile(sourceFigures, 'random_order');
end
destination = fullfile(scriptDirectory, 'paper_results', ...
    config.experimentMode, config.defaultOrdering);
if ~exist(destination, 'dir'), mkdir(destination); end
% Preserve descriptive experimental filenames; README maps them to article panels.
for entry = ["data", "figures"]
    target = fullfile(destination, entry);
    if ~exist(target, 'dir'), mkdir(target); end
    if entry == "data", source = sourceData; else, source = sourceFigures; end
    files = dir(fullfile(source, '*'));
    for index = 1:numel(files)
        if ~files(index).isdir
            copyfile(fullfile(source, files(index).name), target);
        end
    end
    comparison = fullfile(fileparts(source), 'removal_comparison');
    if config.defaultOrdering ~= "random-order"
        comparison = fullfile(source, 'removal_comparison');
    end
    if exist(comparison, 'dir')
        copyfile(comparison, fullfile(target, 'removal_comparison'));
    end
end
copyfile(fullfile(scriptDirectory, 'README.md'), destination);
disp(destination);
