%% Replot existing TSV results only; does not generate matrices or time algorithms.
scriptDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDirectory, 'matlab'));
config = sd_experiment_config();
figureFiles = generate_experiment_figures(config.experimentMode, ...
    config.experimentSwitches, config.defaultOrdering);
disp(figureFiles);
