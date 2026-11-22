function results = run_lambda4_stage_experiment( ...
        mode, family, repetitions, outputFile, ambiguityOrdering)
%RUN_LAMBDA4_STAGE_EXPERIMENT Experiment 1: LAMBDA 4.0 stages only.
%   FAMILY is "random-dimension" or "random-condition". Cached matrices
%   and float ambiguities are read without generating or modifying data.

arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    family (1, 1) string {mustBeMember(family, ...
        ["random-dimension", "random-condition"])} = "random-dimension"
    repetitions (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
    outputFile (1, 1) string = ""
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

paths = setup_sd_experiment_paths();
config = sd_experiment_config(mode, ambiguityOrdering);
if repetitions == 0
    repetitions = config.timingRepetitions;
end
cacheDirectory = paths.cache;
if ambiguityOrdering == "random-order"
    cacheDirectory = fullfile(cacheDirectory, 'random_order');
end
manifestFile = fullfile(cacheDirectory, "manifest_" + mode + ".tsv");
if ~exist(manifestFile, 'file')
    error('run_lambda4_stage_experiment:MissingCache', ...
        'Prepare the %s matrix cache before running experiments.', mode);
end
manifest = readtable(manifestFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
settings = manifest(string(manifest.Family) == family, :);
if family == "random-dimension"
    settings = settings(ismember(settings.Dimension, ...
        config.stageDimensionValues), :);
end
if isempty(settings)
    error('run_lambda4_stage_experiment:MissingFamily', ...
        'No cached settings were found for %s.', family);
end
if strlength(outputFile) == 0
    outputDirectory = fullfile(paths.numericalOutput, mode);
    if ambiguityOrdering == "random-order"
        outputDirectory = fullfile(outputDirectory, 'random_order');
    end
    if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
    if family == "random-dimension"
        fileName = "exp01_lambda4_stage_dimension.tsv";
    else
        fileName = "exp01_lambda4_stage_condition.tsv";
    end
    outputFile = fullfile(outputDirectory, fileName);
else
    outputDirectory = fileparts(outputFile);
    if strlength(outputDirectory) > 0 && ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory);
    end
end

results = table('Size', [height(settings) 8], ...
    'VariableTypes', {'string','double','double','double', ...
                      'double','double','double','double'}, ...
    'VariableNames', {'Family','Dimension','ConditionNumber','NumMatrices', ...
                      'ReductionMeanMs','SearchMeanMs','TotalMeanMs', ...
                      'ReductionFraction'});

for settingIndex = 1:height(settings)
    fprintf('LAMBDA stage experiment: %s setting %d/%d\n', ...
        family, settingIndex, height(settings));
    cacheFile = char(string(settings.File(settingIndex)));
    cached = load(cacheFile, ...
        'matrices', 'floatAmbiguities');
    numMatrices = size(cached.matrices, 3);
    reductionMs = zeros(numMatrices, 1);
    searchMs = zeros(numMatrices, 1);

    % Discard one complete benchmark call at every dimension/condition.
    % MATLAB can continue tiered compilation until a function has returned,
    % which cannot be removed reliably by warm-ups inside the timer wrapper.
    benchmark_lambda4_stages(cached.matrices(:, :, 1), ...
        cached.floatAmbiguities(:, 1), repetitions);
    for sampleIndex = 1:numMatrices
        timing = benchmark_lambda4_stages( ...
            cached.matrices(:, :, sampleIndex), ...
            cached.floatAmbiguities(:, sampleIndex), repetitions);
        reductionMs(sampleIndex) = 1000 * timing.reductionMeanSeconds;
        searchMs(sampleIndex) = 1000 * timing.searchMeanSeconds;
    end
    totalMs = reductionMs + searchMs;
    results(settingIndex, :) = {family, settings.Dimension(settingIndex), ...
        settings.TargetCondition(settingIndex), numMatrices, ...
        mean(reductionMs), mean(searchMs), mean(totalMs), ...
        mean(reductionMs) / mean(totalMs)};
    if strlength(outputFile) > 0
        writetable(results(1:settingIndex, :), outputFile, ...
            'FileType', 'text', 'Delimiter', '\t');
    end
end

writetable(results, outputFile, 'FileType', 'text', 'Delimiter', '\t');

end
