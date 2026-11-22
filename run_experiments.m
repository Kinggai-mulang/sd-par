%% Sequential-decorrelation experiments: user settings
%% Initialization
scriptDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDirectory, 'matlab'));
defaults = sd_experiment_config();
experimentSwitches = defaults.experimentSwitches;
useCachedData = defaults.useCachedData;
experimentMode = defaults.experimentMode;
ambiguityOrdering = defaults.defaultOrdering;
paths = setup_sd_experiment_paths();
config = sd_experiment_config(experimentMode, ambiguityOrdering);
selectedFamilies = strings(0, 1);
if experimentSwitches(1) || experimentSwitches(2) || ...
        experimentSwitches(5) || experimentSwitches(6)
    selectedFamilies(end+1) = "random-dimension";
end
if experimentSwitches(3)
    selectedFamilies(end+1) = "random-condition";
end
if experimentSwitches(4)
    selectedFamilies(end+1) = "geometry-based";
end
if ~any(experimentSwitches)
    error('run_experiments:NoExperiment', ...
        'Enable at least one entry in experimentSwitches.');
end

results = struct();
results.config = config;
if ~isempty(selectedFamilies)
    results.manifest = prepare_matrix_cache( ...
        experimentMode, "", ~useCachedData, ...
        selectedFamilies, ambiguityOrdering);
end
numericalDirectory = fullfile(paths.numericalOutput, experimentMode);
if ambiguityOrdering == "random-order"
    numericalDirectory = fullfile(numericalDirectory, 'random_order');
end
if ~exist(numericalDirectory, 'dir'), mkdir(numericalDirectory); end
if ~any(experimentSwitches(1:4))
    if experimentSwitches(5)
        save(fullfile(numericalDirectory, 'exp05_config.mat'), ...
            'config', 'experimentSwitches', 'useCachedData', 'ambiguityOrdering');
    end
    if experimentSwitches(6)
        save(fullfile(numericalDirectory, 'exp06_random_dimension_par_config.mat'), ...
            'config', 'experimentSwitches', 'useCachedData', 'ambiguityOrdering');
    end
else
    save(fullfile(numericalDirectory, 'experiment_config.mat'), ...
        'config', 'experimentSwitches', 'useCachedData', 'ambiguityOrdering');
end

%% Selected numerical experiments
if experimentSwitches(1)
    results.experiment1 = run_lambda4_stage_experiment( ...
        experimentMode, "random-dimension", config.timingRepetitions, "", ...
        ambiguityOrdering);
end
if experimentSwitches(2)
    results.experiment2 = run_cd_sd_reduction_experiment( ...
        experimentMode, "random-dimension", ...
        config.timingRepetitions, true, "", ambiguityOrdering);
end
if experimentSwitches(3)
    results.experiment3 = run_cd_sd_reduction_experiment( ...
        experimentMode, "random-condition", ...
        config.timingRepetitions, true, "", ambiguityOrdering);
end
if experimentSwitches(4)
    results.experiment4 = run_gb_experiment( ...
        experimentMode, config.timingRepetitions, "", ambiguityOrdering);
end
if experimentSwitches(5)
    results.experiment5 = run_operation_count_experiment( ...
        experimentMode, ambiguityOrdering);
end
if experimentSwitches(6)
    results.experiment6 = run_deletion_chain_experiment( ...
        experimentMode, config.timingRepetitions, "", ...
        ambiguityOrdering, config.numRandomMatrices);
end

if experimentSwitches(7)
    prepare_matrix_cache(experimentMode, "", ~useCachedData, ...
        "random-dimension", "random-order");
    results.experiment7 = run_removal_policy_experiment( ...
        experimentMode, defaults.useCachedResults);
end

%% Figures
results.figureFiles = generate_experiment_figures( ...
    experimentMode, experimentSwitches, ambiguityOrdering);
disp('Selected experiments completed successfully.');
disp(results.figureFiles);
