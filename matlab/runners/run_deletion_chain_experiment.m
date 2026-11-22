function results = run_deletion_chain_experiment( ...
        mode, repetitions, outputFile, ambiguityOrdering, numMatrices, family, cacheRoot, resumeExisting)
%RUN_DELETION_CHAIN_EXPERIMENT Compare complete CD and SD PAR chains.
%   Uses cached random-LDLT or GPS geometry-based inputs. Each chain
%   processes every leading subset n:-1:1 with the same timing code.

arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    repetitions (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
    outputFile (1, 1) string = ""
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
    numMatrices (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
    family (1, 1) string {mustBeMember(family, ["random-dimension", "geometry-based"])} = "random-dimension"
    cacheRoot (1, 1) string = ""
    resumeExisting (1, 1) logical = false
end

paths = setup_sd_experiment_paths();
config = sd_experiment_config(mode, ambiguityOrdering);
if repetitions == 0, repetitions = config.timingRepetitions; end
if numMatrices == 0, numMatrices = config.numRandomMatrices; end

cacheDirectory = paths.cache;
if strlength(cacheRoot) > 0, cacheDirectory = cacheRoot; end
outputDirectory = fullfile(paths.numericalOutput, mode);
if ambiguityOrdering == "random-order"
    cacheDirectory = fullfile(cacheDirectory, 'random_order');
    outputDirectory = fullfile(outputDirectory, 'random_order');
end
manifestFile = fullfile(cacheDirectory, "manifest_" + mode + ".tsv");
if ~exist(manifestFile, 'file')
    error('run_deletion_chain_experiment:MissingCache', ...
        'Generate the requested covariance family with prepare_matrix_cache first.');
end
manifest = readtable(manifestFile, 'FileType', 'text', ...
    'Delimiter', '\t', 'VariableNamingRule', 'preserve');
dimensions = config.dimensionValues;
if family == "geometry-based"
    dimensions = config.gb.numFrequencies * (config.gb.numSatellites - 1);
end
settings = manifest(string(manifest.Family) == family & ...
    ismember(manifest.Dimension, dimensions), :);
settings = sortrows(settings, 'Dimension');
if numMatrices < 1
    error('run_deletion_chain_experiment:InvalidSampleCount', ...
        'At least one covariance matrix is required.');
end
if ~isequal(settings.Dimension(:), dimensions(:))
    error('run_deletion_chain_experiment:IncompleteDimensionCache', ...
        'The %s cache does not contain every configured dimension.', family);
end

if strlength(outputFile) == 0
    if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
    outputFile = fullfile(outputDirectory, ...
        'exp06_random_dimension_full_par.tsv');
    if family == "geometry-based"
        outputFile = fullfile(outputDirectory, 'exp06_gb_full_par.tsv');
    end
else
    explicitDirectory = fileparts(outputFile);
    if strlength(explicitDirectory) > 0 && ~exist(explicitDirectory, 'dir')
        mkdir(explicitDirectory);
    end
end

names = {'Dimension','NumCovarianceMatrices','NumChainLevels', ...
    'CDReductionMeanMs','SDReductionMeanMs', ...
    'CDSearchChainMeanMs','SDSearchChainMeanMs', ...
    'CDMeanMs','SDMeanMs','TimeSpeedup','SearchTimeRatio', ...
    'SolutionMatchRate','MaxObjectiveDifference'};
results = array2table(NaN(height(settings), numel(names)), ...
    'VariableNames', names);
firstSetting = 1;
if resumeExisting && isfile(outputFile)
    completed = readtable(outputFile, 'FileType', 'text', 'Delimiter', '\t');
    assert(isequal(completed.Properties.VariableNames, names));
    assert(height(completed) <= height(settings));
    assert(isequal(completed.Dimension, settings.Dimension(1:height(completed))));
    assert(all(completed.NumCovarianceMatrices == numMatrices));
    assert(all(completed.NumChainLevels == completed.Dimension));
    assert(all(completed.SolutionMatchRate == 1));
    results(1:height(completed), :) = completed;
    firstSetting = height(completed) + 1;
    fprintf('Resume: retaining %d completed settings.\n', height(completed));
end

for settingIndex = firstSetting:height(settings)
    dimension = settings.Dimension(settingIndex);
    fprintf('Experiment 6 (%s PAR): setting %d/%d, n=%d\n', ...
        family, settingIndex, height(settings), dimension);
    cached = load(char(string(settings.File(settingIndex))), ...
        'matrices', 'floatAmbiguities');
    if size(cached.matrices, 1) ~= dimension || ...
            size(cached.matrices, 2) ~= dimension || ...
            size(cached.matrices, 3) < numMatrices || ...
            size(cached.floatAmbiguities, 2) < numMatrices
        error('run_deletion_chain_experiment:InsufficientCachedSamples', ...
            'Dimension %d cache has fewer than %d paired samples.', ...
            dimension, numMatrices);
    end

    matrices = cached.matrices(:, :, 1:numMatrices);
    floatAmbiguities = cached.floatAmbiguities(:, 1:numMatrices);
    for warmupIndex = 1:min(2, repetitions)
        runCompleteChain(matrices(:, :, 1), ...
            floatAmbiguities(:, 1), false);
        runCompleteChain(matrices(:, :, 1), ...
            floatAmbiguities(:, 1), true);
    end

    cdTotal = zeros(numMatrices, 1);
    sdTotal = zeros(numMatrices, 1);
    cdReduction = zeros(numMatrices, 1);
    sdReduction = zeros(numMatrices, 1);
    cdSearch = zeros(numMatrices, 1);
    sdSearch = zeros(numMatrices, 1);
    solutionMatches = zeros(numMatrices, 1);
    objectiveDifferences = zeros(numMatrices, 1);
    firstSample = 1;
    checkpointFile = outputFile + ".n" + compose("%03d", dimension) + ".checkpoint.mat";
    if resumeExisting && isfile(checkpointFile)
        saved = load(checkpointFile);
        assert(saved.checkpoint.dimension == dimension);
        assert(saved.checkpoint.numMatrices == numMatrices);
        assert(saved.checkpoint.family == family);
        assert(saved.checkpoint.ambiguityOrdering == ambiguityOrdering);
        assert(saved.checkpoint.sourceFile == string(settings.File(settingIndex)));
        assert(saved.checkpoint.completedSamples >= 0 && ...
            saved.checkpoint.completedSamples <= numMatrices && ...
            saved.checkpoint.completedSamples == fix(saved.checkpoint.completedSamples));
        assert(numel(saved.cdTotal) == numMatrices && numel(saved.sdTotal) == numMatrices);
        cdTotal = saved.cdTotal; sdTotal = saved.sdTotal;
        cdReduction = saved.cdReduction; sdReduction = saved.sdReduction;
        cdSearch = saved.cdSearch; sdSearch = saved.sdSearch;
        solutionMatches = saved.solutionMatches;
        objectiveDifferences = saved.objectiveDifferences;
        firstSample = saved.checkpoint.completedSamples + 1;
        fprintf('  Resume samples: retaining %d paired matrices.\n', firstSample - 1);
    end
    progressStep = max(1, floor(numMatrices / 20));
    for sampleIndex = firstSample:numMatrices
        Qa = matrices(:, :, sampleIndex);
        aFloat = floatAmbiguities(:, sampleIndex);
        if mod(settingIndex + sampleIndex, 2) == 0
            startTime = tic;
            cd = runCompleteChain(Qa, aFloat, false);
            cdTotal(sampleIndex) = toc(startTime);
            startTime = tic;
            sd = runCompleteChain(Qa, aFloat, true);
            sdTotal(sampleIndex) = toc(startTime);
        else
            startTime = tic;
            sd = runCompleteChain(Qa, aFloat, true);
            sdTotal(sampleIndex) = toc(startTime);
            startTime = tic;
            cd = runCompleteChain(Qa, aFloat, false);
            cdTotal(sampleIndex) = toc(startTime);
        end
        cdReduction(sampleIndex) = cd.reductionSeconds;
        sdReduction(sampleIndex) = sd.reductionSeconds;
        cdSearch(sampleIndex) = cd.searchSeconds;
        sdSearch(sampleIndex) = sd.searchSeconds;
        solutionMatches(sampleIndex) = mean(all( ...
            cd.fixedAmbiguities == sd.fixedAmbiguities, 1));
        objectiveDifferences(sampleIndex) = max(abs( ...
            cd.objectives - sd.objectives));
        if mod(sampleIndex, progressStep) == 0 || sampleIndex == numMatrices
            if resumeExisting
                checkpoint = struct('dimension', dimension, 'numMatrices', numMatrices, ...
                    'family', family, 'ambiguityOrdering', ambiguityOrdering, ...
                    'sourceFile', string(settings.File(settingIndex)), ...
                    'completedSamples', sampleIndex);
                % Save outside both methods' timing regions, then replace the checkpoint.
                temporaryCheckpoint = checkpointFile + ".tmp.mat";
                save(temporaryCheckpoint, 'checkpoint', 'cdTotal', 'sdTotal', ...
                    'cdReduction', 'sdReduction', 'cdSearch', 'sdSearch', ...
                    'solutionMatches', 'objectiveDifferences');
                movefile(temporaryCheckpoint, checkpointFile, 'f');
            end
            fprintf('  Complete chains: %d/%d paired matrices\n', ...
                sampleIndex, numMatrices);
        end
    end

    cdMeanMs = 1000 * mean(cdTotal);
    sdMeanMs = 1000 * mean(sdTotal);
    cdSearchMeanMs = 1000 * mean(cdSearch);
    sdSearchMeanMs = 1000 * mean(sdSearch);
    results(settingIndex, :) = {dimension, numMatrices, dimension, ...
        1000 * mean(cdReduction), 1000 * mean(sdReduction), ...
        cdSearchMeanMs, sdSearchMeanMs, cdMeanMs, sdMeanMs, ...
        cdMeanMs / sdMeanMs, cdSearchMeanMs / sdSearchMeanMs, ...
        mean(solutionMatches), max(objectiveDifferences)};
    writetable(results(1:settingIndex, :), outputFile, ...
        'FileType', 'text', 'Delimiter', '\t');
end
end

function result = runCompleteChain(Qa, aFloat, useSequential)
n = size(Qa, 1);
result.reductionSeconds = 0;
result.searchSeconds = 0;
result.fixedAmbiguities = zeros(n, n);
result.objectives = zeros(1, n);
if useSequential
    parentCovariance = Qa;
end
for dimension = n:-1:1
    startTime = tic;
    if useSequential
        if dimension == n
            reduction = lambda4_reduce_covariance(Qa);
        else
            reduction = sd_lambda4_reduce_child( ...
                parentCovariance, reduction, dimension + 1);
            parentCovariance = reduction.Qsubset;
        end
    else
        reduction = lambda4_reduce_covariance( ...
            Qa(1:dimension, 1:dimension));
    end
    result.reductionSeconds = result.reductionSeconds + toc(startTime);

    zFloat = reduction.Z * aFloat(1:dimension);
    startTime = tic;
    [zFixed, objective] = estimatorILS( ...
        zFloat, reduction.L, reduction.d, 2);
    result.searchSeconds = result.searchSeconds + toc(startTime);
    aFixed = reduction.inverseZ * zFixed(:, 1);
    if any(~isfinite(aFixed))
        error('run_deletion_chain_experiment:InvalidSearchResult', ...
            'ILS returned a non-finite back-transformed solution.');
    end
    result.fixedAmbiguities(1:dimension, dimension) = aFixed;
    result.objectives(dimension) = objective(1);
end
end
