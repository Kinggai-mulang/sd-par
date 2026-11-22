function results = run_removal_policy_experiment(mode, useCachedResults)
%RUN_REMOVAL_POLICY_EXPERIMENT Paired last/random, max and min variance removal.
%   Reuses Experiment 2 random-order inputs without changing their cache.
arguments
    mode (1,1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    useCachedResults (1,1) logical = true
end
paths = setup_sd_experiment_paths();
config = sd_experiment_config(mode, "random-order");
manifest = read_cached_experiment_manifest( ...
    paths.cache, mode, "random-order", true);
settings = manifest(string(manifest.Family) == "random-dimension" & ...
    ismember(manifest.Dimension, config.dimensionValues), :);
settings = sortrows(settings, 'Dimension');
if ~isequal(settings.Dimension(:)', config.dimensionValues)
    error('run_removal_policy_experiment:IncompleteCache', ...
        'Prepare the Experiment 2 random-order cache for every configured dimension.');
end
outputDirectory = fullfile(paths.numericalOutput, mode, 'removal_comparison');
if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
outputFile = fullfile(outputDirectory, 'removal_policy_dimension.tsv');
protocol = 'paired-removal-stable-survivors-v1';
policies = ["random-removal", "maximum-variance", "minimum-variance"];
save(fullfile(outputDirectory, 'removal_comparison_config.mat'), ...
    'config', 'settings', 'protocol', 'policies');
names = {'Dimension','NumMatrices','TimingRepetitions', ...
    'RandomCDMeanMs','RandomSDMeanMs','RandomSpeedup', ...
    'MaximumVarianceCDMeanMs','MaximumVarianceSDMeanMs','MaximumVarianceSpeedup', ...
    'MinimumVarianceCDMeanMs','MinimumVarianceSDMeanMs','MinimumVarianceSpeedup'};
results = array2table(NaN(height(settings), numel(names)), 'VariableNames', names);

for settingIndex = 1:height(settings)
    n = settings.Dimension(settingIndex);
    fprintf('Removal comparison: setting %d/%d, n=%d\n', ...
        settingIndex, height(settings), n);
    cached = load(char(string(settings.File(settingIndex))), 'matrices','metadata');
    matrices = cached.matrices;
    sourceMetadata = cached.metadata;
    numMatrices = size(matrices, 3);
    repetitions = config.timingRepetitions;
    if numMatrices ~= config.numRandomMatrices || ...
            string(sourceMetadata.ordering) ~= "random-order"
        error('run_removal_policy_experiment:CacheMismatch', ...
            'Input count or ordering differs from Experiment 2 configuration.');
    end
    checkpointFile = fullfile(outputDirectory, sprintf('n%03d_timings.mat', n));
    reused = false;
    if useCachedResults && exist(checkpointFile, 'file')
        checkpoint = load(checkpointFile);
        reused = isequaln(checkpoint.sourceMetadata, sourceMetadata) ...
            && checkpoint.repetitions == repetitions ...
            && checkpoint.numMatrices == numMatrices ...
            && string(checkpoint.protocol) == protocol;
        if reused
            cdMs = checkpoint.cdMs;
            sdMs = checkpoint.sdMs;
        end
    end
    if ~reused
        cdMs = zeros(numMatrices, 3);
        sdMs = zeros(numMatrices, 3);
        removedIndices = zeros(numMatrices, 3);
        % Discard full benchmark calls to exercise all three paths before timing.
        first = matrices(:,:,1);
        [~, maxIndex] = max(diag(first));
        [~, minIndex] = min(diag(first));
        for removeIndex = [n maxIndex minIndex]
            order = [1:removeIndex-1, removeIndex+1:n, removeIndex];
            benchmark_cd_sd_reduction(first(order,order), repetitions);
        end
        progressStep = max(1, floor(numMatrices/10));
        for sampleIndex = 1:numMatrices
            Qa = matrices(:,:,sampleIndex);
            [~, maxIndex] = max(diag(Qa));
            [~, minIndex] = min(diag(Qa));
            removedIndices(sampleIndex,:) = [n maxIndex minIndex];
            % Rotate policy execution order across samples to limit timing bias.
            policyOrder = mod((0:2)+sampleIndex-1, 3)+1;
            for policyIndex = policyOrder
                removeIndex = removedIndices(sampleIndex,policyIndex);
                order = [1:removeIndex-1, removeIndex+1:n, removeIndex];
                try
                    timing = benchmark_cd_sd_reduction(Qa(order,order), repetitions);
                catch cause
                    error('run_removal_policy_experiment:ReductionFailure', ...
                        'n=%d, sample=%d, policy=%s: %s', ...
                        n, sampleIndex, policies(policyIndex), cause.message);
                end
                cdMs(sampleIndex,policyIndex) = 1000*timing.cdMeanSeconds;
                sdMs(sampleIndex,policyIndex) = 1000*timing.sdMeanSeconds;
            end
            if mod(sampleIndex, progressStep) == 0
                fprintf('  paired matrices %d/%d\n', sampleIndex, numMatrices);
            end
        end
        save(checkpointFile, 'cdMs','sdMs','removedIndices','sourceMetadata', ...
            'repetitions','numMatrices','protocol','policies');
    end
    cdMean = mean(cdMs,1);
    sdMean = mean(sdMs,1);
    ratios = cdMean./sdMean;
    results(settingIndex,:) = num2cell([n numMatrices repetitions, ...
        cdMean(1) sdMean(1) ratios(1), cdMean(2) sdMean(2) ratios(2), ...
        cdMean(3) sdMean(3) ratios(3)]);
    writetable(results(1:settingIndex,:), outputFile, ...
        'FileType','text','Delimiter','\t');
end
end
