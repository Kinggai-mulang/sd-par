function manifest = prepare_matrix_cache( ...
        mode, cacheDirectory, overwrite, families, ambiguityOrdering)
%PREPARE_MATRIX_CACHE Generate, validate, and save reusable input matrices.
%   One MAT file is written per parameter setting. Existing files are
%   validated and reused unless OVERWRITE is true.

arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    cacheDirectory (1, 1) string = ""
    overwrite (1, 1) logical = false
    families (1, :) string {mustBeMember(families, ...
        ["random-dimension", "random-condition", "geometry-based"])} = ...
        ["random-dimension", "random-condition", "geometry-based"]
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

paths = setup_sd_experiment_paths();
config = sd_experiment_config(mode, ambiguityOrdering);
if strlength(cacheDirectory) == 0
    cacheDirectory = string(paths.cache);
    if ambiguityOrdering == "random-order"
        cacheDirectory = fullfile(cacheDirectory, 'random_order');
    end
end

dimensionDirectory = fullfile(cacheDirectory, 'random_dimension');
conditionDirectory = fullfile(cacheDirectory, 'random_condition');
gbDirectory = fullfile(cacheDirectory, 'geometry_based');
if ~exist(dimensionDirectory, 'dir'), mkdir(dimensionDirectory); end
if ~exist(conditionDirectory, 'dir'), mkdir(conditionDirectory); end
if ~exist(gbDirectory, 'dir'), mkdir(gbDirectory); end

manifestFile = fullfile(cacheDirectory, "manifest_" + mode + ".tsv");
manifest = table('Size', [0 7], ...
    'VariableTypes', {'string','string','double','double','double','double','string'}, ...
    'VariableNames', {'Family','File','Dimension','TargetCondition', ...
                      'ActualConditionMin','ActualConditionMax','Status'});
if exist(manifestFile, 'file')
    previous = readtable(manifestFile, 'FileType', 'text', ...
        'Delimiter', '\t', 'VariableNamingRule', 'preserve');
    previous.Family = string(previous.Family);
    previous.File = string(previous.File);
    previous.Status = string(previous.Status);
    manifest = previous(~ismember(string(previous.Family), families), :);
end

if ismember("random-dimension", families)
dimensionValues = unique([config.stageDimensionValues, config.dimensionValues]);
for settingIndex = 1:numel(dimensionValues)
    n = dimensionValues(settingIndex);
    kappa = NaN;
    file = fullfile(dimensionDirectory, sprintf('n%03d_ldlt_%s.mat', n, mode));
    seedBase = config.baseSeed + 100000 * settingIndex;
    [actualMin, actualMax, status] = prepareRandomFile( ...
        file, n, kappa, seedBase, config.numRandomMatrices, ...
        "random-ldlt", overwrite, ambiguityOrdering);
    manifest(end+1, :) = {'random-dimension', string(file), n, kappa, ...
        actualMin, actualMax, status}; %#ok<AGROW>
end
end

if ismember("random-condition", families)
for settingIndex = 1:numel(config.conditionExponents)
    exponent = config.conditionExponents(settingIndex);
    n = config.conditionDimension;
    kappa = 2^exponent;
    file = fullfile(conditionDirectory, sprintf('n%03d_k2p%02d_%s.mat', ...
        n, exponent, mode));
    seedBase = config.baseSeed + 2000000 + 100000 * settingIndex;
    [actualMin, actualMax, status] = prepareRandomFile( ...
        file, n, kappa, seedBase, config.numRandomMatrices, ...
        "controlled-spd", overwrite, ambiguityOrdering);
    manifest(end+1, :) = {'random-condition', string(file), n, kappa, ...
        actualMin, actualMax, status}; %#ok<AGROW>
end
end

if ismember("geometry-based", families)
for settingIndex = 1:numel(config.gb.numSatellites)
    numSatellites = config.gb.numSatellites(settingIndex);
    n = config.gb.numFrequencies * (numSatellites - 1);
    file = fullfile(gbDirectory, sprintf( ...
        'gps_gb_r2_f3_s%02d_n%03d_%s.mat', numSatellites, n, mode));
    seedBase = config.baseSeed + 3000000 + 100000 * settingIndex;
    needsGeneration = overwrite || ~exist(file, 'file');
    if ~needsGeneration
        cached = load(file);
        needsGeneration = ~isfield(cached, 'metadata') ...
            || ~isfield(cached.metadata, 'generator') ...
            || string(cached.metadata.generator) ~= "gps-geometry-based" ...
            || ~isfield(cached.metadata, 'ordering') ...
            || string(cached.metadata.ordering) ~= ambiguityOrdering ...
            || ~isfield(cached.metadata, 'seedBase') ...
            || cached.metadata.seedBase ~= seedBase ...
            || ~isfield(cached.metadata, 'settings') ...
            || ~isequaln(cached.metadata.settings, config.gb) ...
            || ~isfield(cached, 'matrices') ...
            || size(cached.matrices, 1) ~= n ...
            || size(cached.matrices, 2) ~= n ...
            || size(cached.matrices, 3) ~= config.numRandomMatrices ...
            || ~isfield(cached, 'floatAmbiguities') ...
            || ~isequal(size(cached.floatAmbiguities), ...
                [n config.numRandomMatrices]);
    end
    if needsGeneration
        matrices = zeros(n, n, config.numRandomMatrices);
        floatAmbiguities = zeros(n, config.numRandomMatrices);
        actualConditions = zeros(config.numRandomMatrices, 1);
        geometrySeeds = zeros(config.numRandomMatrices, 1);
        permutations = zeros(n, config.numRandomMatrices);
        for sampleIndex = 1:config.numRandomMatrices
            candidateSeed = seedBase + sampleIndex - 1;
            for attempt = 1:100
                try
                    [Qa, generated] = generate_gb_covariance( ...
                        config.gb, numSatellites, candidateSeed, ...
                        ambiguityOrdering);
                    break
                catch cause
                    if ~ismember(string(cause.identifier), [ ...
                            "validate_covariance_matrix:LambdaDecompositionFailure", ...
                            "generate_gb_covariance:SingularGeometry"])
                        rethrow(cause)
                    end
                    if attempt == 100
                        rethrow(cause)
                    end
                    candidateSeed = candidateSeed + 10000000;
                end
            end
            matrices(:, :, sampleIndex) = Qa;
            actualConditions(sampleIndex) = generated.actualCondition;
            geometrySeeds(sampleIndex) = candidateSeed;
            permutations(:, sampleIndex) = generated.permutation(:);
            stream = RandStream('mt19937ar', ...
                'Seed', seedBase + 50000000 + sampleIndex - 1);
            floatAmbiguities(:, sampleIndex) = ...
                chol(Qa, 'lower') * randn(stream, n, 1);
        end
        metadata.generator = 'gps-geometry-based';
        metadata.settings = config.gb;
        metadata.numSatellites = numSatellites;
        metadata.numMatrices = config.numRandomMatrices;
        metadata.seedBase = seedBase;
        metadata.geometrySeeds = geometrySeeds;
        metadata.permutations = permutations;
        metadata.actualConditions = actualConditions;
        metadata.ordering = char(ambiguityOrdering);
        metadata.removal = 'remove last after ordering';
        save(file, 'matrices', 'floatAmbiguities', 'metadata', '-v7.3');
        status = "generated";
    else
        validate_covariance_matrix(cached.matrices(:, :, 1), ...
            n, [], ambiguityOrdering);
        validate_covariance_matrix(cached.matrices(:, :, end), ...
            n, [], ambiguityOrdering);
        actualConditions = cached.metadata.actualConditions;
        status = "reused";
    end
    manifest(end+1, :) = {'geometry-based', string(file), n, NaN, ...
        min(actualConditions), max(actualConditions), status}; %#ok<AGROW>
end
end

manifest = sortrows(manifest, {'Family', 'Dimension'});
writetable(manifest, manifestFile, ...
    'FileType', 'text', 'Delimiter', '\t');

end

function [actualMin, actualMax, status] = prepareRandomFile( ...
        file, dimension, targetCondition, seedBase, numMatrices, ...
        generatorKind, overwrite, ambiguityOrdering)

needsGeneration = overwrite || ~exist(file, 'file');
if ~needsGeneration
    cached = load(file);
    cachedOrdering = "";
    if isfield(cached, 'metadata') && isfield(cached.metadata, 'ordering')
        cachedOrdering = string(cached.metadata.ordering);
    end
    legacyVarianceOrder = ambiguityOrdering == "variance-ordered" ...
        && cachedOrdering == "ascending marginal variance; remove last";
    needsGeneration = ~isfield(cached, 'metadata') ...
        || ~isfield(cached.metadata, 'generator') ...
        || string(cached.metadata.generator) ~= generatorKind ...
        || ~isfield(cached.metadata, 'ordering') ...
        || ~(cachedOrdering == ambiguityOrdering || legacyVarianceOrder) ...
        || ~isfield(cached.metadata, 'seedBase') ...
        || cached.metadata.seedBase ~= seedBase ...
        || ~isfield(cached.metadata, 'targetCondition') ...
        || ~isequaln(cached.metadata.targetCondition, targetCondition) ...
        || ~isfield(cached, 'matrices') ...
        || size(cached.matrices, 1) ~= dimension ...
        || size(cached.matrices, 2) ~= dimension ...
        || size(cached.matrices, 3) ~= numMatrices;
end

if ~needsGeneration
    if size(cached.matrices, 3) ~= numMatrices
        error('prepare_matrix_cache:WrongCachedCount', ...
            'Cached file %s contains %d matrices instead of %d.', ...
            file, size(cached.matrices, 3), numMatrices);
    end
    validate_covariance_matrix(cached.matrices(:, :, 1), ...
        dimension, targetCondition, ambiguityOrdering);
    validate_covariance_matrix(cached.matrices(:, :, end), ...
        dimension, targetCondition, ambiguityOrdering);
    actualConditions = cached.metadata.actualConditions;
    if isfield(cached, 'floatAmbiguities') && ...
            isequal(size(cached.floatAmbiguities), [dimension numMatrices])
        status = "reused";
    else
        matrices = cached.matrices;
        floatAmbiguities = generateFloatAmbiguities( ...
            matrices, seedBase + 50000000);
        metadata = cached.metadata;
        save(file, 'matrices', 'floatAmbiguities', 'metadata');
        status = "upgraded";
    end
else
    matrices = zeros(dimension, dimension, numMatrices);
    actualConditions = zeros(numMatrices, 1);
    matrixSeeds = zeros(numMatrices, 1);
    for sampleIndex = 1:numMatrices
        candidateSeed = seedBase + sampleIndex - 1;
        for attempt = 1:100
            try
                if generatorKind == "random-ldlt"
                    [matrices(:, :, sampleIndex), generated] = ...
                        generate_random_ldlt( ...
                        dimension, candidateSeed, ambiguityOrdering);
                else
                    [matrices(:, :, sampleIndex), generated] = ...
                        generate_controlled_spd( ...
                        dimension, targetCondition, candidateSeed, ...
                        ambiguityOrdering);
                end
                break
            catch cause
                if string(cause.identifier) ~= ...
                        "validate_covariance_matrix:LambdaDecompositionFailure"
                    rethrow(cause)
                end
                if attempt == 100
                    error('prepare_matrix_cache:NoStableMatrix', ...
                        ['Could not generate a LAMBDA-compatible matrix ', ...
                         'after 100 deterministic attempts.']);
                end
                candidateSeed = candidateSeed + 10000000;
            end
        end
        actualConditions(sampleIndex) = generated.actualCondition;
        matrixSeeds(sampleIndex) = generated.seed;
    end
    floatAmbiguities = generateFloatAmbiguities( ...
        matrices, seedBase + 50000000);
    metadata.generator = char(generatorKind);
    metadata.dimension = dimension;
    metadata.targetCondition = targetCondition;
    metadata.seedBase = seedBase;
    metadata.numMatrices = numMatrices;
    metadata.actualConditions = actualConditions;
    metadata.matrixSeeds = matrixSeeds;
    metadata.ordering = char(ambiguityOrdering);
    metadata.removal = 'remove last after ordering';
    save(file, 'matrices', 'floatAmbiguities', 'metadata');
    status = "generated";
end

actualMin = min(actualConditions);
actualMax = max(actualConditions);

end

function floatAmbiguities = generateFloatAmbiguities(matrices, seedBase)
dimension = size(matrices, 1);
numMatrices = size(matrices, 3);
floatAmbiguities = zeros(dimension, numMatrices);
for sampleIndex = 1:numMatrices
    stream = RandStream('mt19937ar', 'Seed', seedBase + sampleIndex - 1);
    lowerCholesky = chol(matrices(:, :, sampleIndex), 'lower');
    floatAmbiguities(:, sampleIndex) = lowerCholesky ...
        * randn(stream, dimension, 1);
end
end
