function results = run_cd_sd_reduction_experiment( ...
        mode, family, repetitions, includeFlops, outputFile, ambiguityOrdering)
%RUN_CD_SD_REDUCTION_EXPERIMENT Experiment 2: CD versus SD reduction.
%   Timing always uses uninstrumented implementations. FLOPs are evaluated
%   in a separate call.

arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    family (1, 1) string {mustBeMember(family, ...
        ["random-dimension", "random-condition"])} = "random-dimension"
    repetitions (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
    includeFlops (1, 1) logical = true
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
    error('run_cd_sd_reduction_experiment:MissingCache', ...
        'Prepare the %s matrix cache before running experiments.', mode);
end
manifest = readtable(manifestFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
settings = manifest(string(manifest.Family) == family, :);
if family == "random-dimension"
    settings = settings(ismember(settings.Dimension, ...
        config.dimensionValues), :);
end
if isempty(settings)
    error('run_cd_sd_reduction_experiment:MissingFamily', ...
        'No cached settings were found for %s.', family);
end
if strlength(outputFile) == 0
    outputDirectory = fullfile(paths.numericalOutput, mode);
    if ambiguityOrdering == "random-order"
        outputDirectory = fullfile(outputDirectory, 'random_order');
    end
    if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
    switch family
        case "random-dimension"
            fileName = "exp02_cd_sd_dimension.tsv";
        case "random-condition"
            fileName = "exp03_cd_sd_condition.tsv";
    end
    outputFile = fullfile(outputDirectory, fileName);
else
    outputDirectory = fileparts(outputFile);
    if strlength(outputDirectory) > 0 && ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory);
    end
end

variableNames = {'Family','Dimension','ConditionNumber','NumMatrices', ...
    'CDMeanMs','SDMeanMs','TimeSpeedup', ...
    'CDMeanFlops','SDMeanFlops','FlopSpeedup'};
variableTypes = [{'string'}, repmat({'double'}, 1, numel(variableNames)-1)];
results = table('Size', [height(settings) numel(variableNames)], ...
    'VariableTypes', variableTypes, 'VariableNames', variableNames);

for settingIndex = 1:height(settings)
    fprintf('CD/SD reduction experiment: %s setting %d/%d\n', ...
        family, settingIndex, height(settings));
    cacheFile = char(string(settings.File(settingIndex)));
    cached = load(cacheFile);
    if isfield(cached, 'matrices')
        matrices = cached.matrices;
    else
        matrices = reshape(cached.Qa, size(cached.Qa, 1), size(cached.Qa, 2), 1);
    end
    numMatrices = size(matrices, 3);

    cdMs = zeros(numMatrices, 1);
    sdMs = zeros(numMatrices, 1);
    cdFlops = NaN(numMatrices, 1);
    sdFlops = NaN(numMatrices, 1);

    % Discard one complete benchmark call for this parameter setting. This
    % removes function-level tiered compilation from the first real sample.
    benchmark_cd_sd_reduction(matrices(:, :, 1), repetitions);
    for sampleIndex = 1:numMatrices
        Qa = matrices(:, :, sampleIndex);
        timing = benchmark_cd_sd_reduction(Qa, repetitions);
        cdMs(sampleIndex) = 1000 * timing.cdMeanSeconds;
        sdMs(sampleIndex) = 1000 * timing.sdMeanSeconds;

        if includeFlops
            counts = count_cd_sd_reduction_flops(Qa);
            cdFlops(sampleIndex) = counts.cd;
            sdFlops(sampleIndex) = counts.sd;
        end
    end

    results(settingIndex, :) = {family, settings.Dimension(settingIndex), ...
        settings.TargetCondition(settingIndex), numMatrices, ...
        mean(cdMs), mean(sdMs), mean(cdMs) / mean(sdMs), ...
        mean(cdFlops), mean(sdFlops), mean(cdFlops) / mean(sdFlops)};
    if strlength(outputFile) > 0
        writetable(results(1:settingIndex, :), outputFile, ...
            'FileType', 'text', 'Delimiter', '\t');
    end
end

writetable(results, outputFile, 'FileType', 'text', 'Delimiter', '\t');

end
