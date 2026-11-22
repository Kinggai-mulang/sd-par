function results = run_gb_experiment( ...
        mode, repetitions, outputFile, ambiguityOrdering)
%RUN_GB_EXPERIMENT Compare CD/SD on independent GB covariance matrices.

arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    repetitions (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
    outputFile (1, 1) string = ""
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

paths = setup_sd_experiment_paths();
config = sd_experiment_config(mode, ambiguityOrdering);
if repetitions == 0, repetitions = config.timingRepetitions; end
cacheDirectory = paths.cache;
if ambiguityOrdering == "random-order"
    cacheDirectory = fullfile(cacheDirectory, 'random_order');
end
manifest = readtable(fullfile(cacheDirectory, "manifest_" + mode + ".tsv"), ...
    'FileType', 'text', 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
settings = manifest(string(manifest.Family) == "geometry-based", :);
if height(settings) ~= numel(config.gb.numSatellites)
    error('run_gb_experiment:MissingCache', ...
        'Prepare the complete geometry-based cache before Experiment 4.');
end
if strlength(outputFile) == 0
    outputDirectory = fullfile(paths.numericalOutput, mode);
    if ambiguityOrdering == "random-order"
        outputDirectory = fullfile(outputDirectory, 'random_order');
    end
    if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
    outputFile = fullfile(outputDirectory, 'exp04_gb.tsv');
end

names = {'Dimension','NumSatellites','NumCovarianceMatrices', ...
    'NumSearchSamples','CDMeanMs','SDMeanMs','TimeSpeedup', ...
    'CDSearchMeanMs','SDSearchMeanMs','SearchTimeRatio', ...
    'SolutionMatchRate','MaxObjectiveDifference'};
results = array2table(NaN(height(settings), numel(names)), ...
    'VariableNames', names);

for settingIndex = 1:height(settings)
    fprintf('Experiment 4 (GB): setting %d/%d\n', ...
        settingIndex, height(settings));
    cached = load(char(string(settings.File(settingIndex))), ...
        'matrices', 'floatAmbiguities', 'metadata');
    matrices = cached.matrices;
    floatAmbiguities = cached.floatAmbiguities;
    n = size(matrices, 1);
    numMatrices = size(matrices, 3);
    if numMatrices ~= config.numRandomMatrices || ...
            ~isequal(size(floatAmbiguities), [n numMatrices]) || ...
            string(cached.metadata.ordering) ~= ambiguityOrdering
        error('run_gb_experiment:InvalidCache', ...
            'GB sample count or ordering does not match the configuration.');
    end

    cdReductionSeconds = zeros(numMatrices, 1);
    sdReductionSeconds = zeros(numMatrices, 1);
    cdFloat = cell(numMatrices, 1);
    sdFloat = cell(numMatrices, 1);
    cdL = cell(numMatrices, 1);
    sdL = cell(numMatrices, 1);
    cdD = cell(numMatrices, 1);
    sdD = cell(numMatrices, 1);
    cdInverseZ = cell(numMatrices, 1);
    sdInverseZ = cell(numMatrices, 1);
    for sampleIndex = 1:numMatrices
        Qa = matrices(:, :, sampleIndex);
        timing = benchmark_cd_sd_reduction(Qa, repetitions);
        cdReductionSeconds(sampleIndex) = timing.cdMeanSeconds;
        sdReductionSeconds(sampleIndex) = timing.sdMeanSeconds;
        Qchild = Qa(1:n-1, 1:n-1);
        parent = lambda4_reduce_covariance(Qa);
        cd = lambda4_reduce_covariance(Qchild);
        sd = sd_lambda4_reduce_child(Qa, parent, n);
        childFloat = floatAmbiguities(1:n-1, sampleIndex);
        cdFloat{sampleIndex} = cd.Z*childFloat;
        sdFloat{sampleIndex} = sd.Z*childFloat;
        cdL{sampleIndex} = cd.L;
        sdL{sampleIndex} = sd.L;
        cdD{sampleIndex} = cd.d;
        sdD{sampleIndex} = sd.d;
        cdInverseZ{sampleIndex} = cd.inverseZ;
        sdInverseZ{sampleIndex} = sd.inverseZ;
    end

    % Warm each ILS path, then time whole batches of different covariances.
    for warmupIndex = 1:2
        estimatorILS(cdFloat{1}, cdL{1}, cdD{1}, 2);
        estimatorILS(sdFloat{1}, sdL{1}, sdD{1}, 2);
    end
    cdBatchSeconds = zeros(config.gbSearchBatchRepetitions, 1);
    sdBatchSeconds = zeros(config.gbSearchBatchRepetitions, 1);
    for batchIndex = 1:config.gbSearchBatchRepetitions
        if mod(settingIndex+batchIndex, 2) == 0
            startTime = tic;
            runSearchBatch(cdFloat, cdL, cdD);
            cdBatchSeconds(batchIndex) = toc(startTime);
            startTime = tic;
            runSearchBatch(sdFloat, sdL, sdD);
            sdBatchSeconds(batchIndex) = toc(startTime);
        else
            startTime = tic;
            runSearchBatch(sdFloat, sdL, sdD);
            sdBatchSeconds(batchIndex) = toc(startTime);
            startTime = tic;
            runSearchBatch(cdFloat, cdL, cdD);
            cdBatchSeconds(batchIndex) = toc(startTime);
        end
    end
    matches = false(numMatrices, 1);
    objectiveDifference = zeros(numMatrices, 1);
    for sampleIndex = 1:numMatrices
        [cdZFix, cdNorm] = estimatorILS( ...
            cdFloat{sampleIndex}, cdL{sampleIndex}, cdD{sampleIndex}, 2);
        [sdZFix, sdNorm] = estimatorILS( ...
            sdFloat{sampleIndex}, sdL{sampleIndex}, sdD{sampleIndex}, 2);
        matches(sampleIndex) = isequal( ...
            cdInverseZ{sampleIndex}*cdZFix(:, 1), ...
            sdInverseZ{sampleIndex}*sdZFix(:, 1));
        objectiveDifference(sampleIndex) = abs(cdNorm(1)-sdNorm(1));
    end

    cdSearchSeconds = median(cdBatchSeconds)/numMatrices;
    sdSearchSeconds = median(sdBatchSeconds)/numMatrices;
    cdReductionMean = mean(cdReductionSeconds);
    sdReductionMean = mean(sdReductionSeconds);
    results(settingIndex, :) = {n, cached.metadata.numSatellites, ...
        numMatrices, numMatrices, 1000*cdReductionMean, ...
        1000*sdReductionMean, cdReductionMean/sdReductionMean, ...
        1000*cdSearchSeconds, 1000*sdSearchSeconds, ...
        cdSearchSeconds/sdSearchSeconds, mean(matches), ...
        max(objectiveDifference)};
    writetable(results(1:settingIndex, :), outputFile, ...
        'FileType', 'text', 'Delimiter', '\t');
end
end

function runSearchBatch(floatAmbiguities, L, d)
for sampleIndex = 1:numel(floatAmbiguities)
    estimatorILS(floatAmbiguities{sampleIndex}, ...
        L{sampleIndex}, d{sampleIndex}, 2);
end
end
