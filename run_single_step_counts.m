function summary = run_single_step_counts(family, mode)
%RUN_SINGLE_STEP_COUNTS Count reduction work on generated covariance inputs.
arguments
    family (1, 1) string {mustBeMember(family, ...
        ["random-dimension", "random-condition", "geometry-based"])} = "geometry-based"
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
end

directory = fileparts(mfilename('fullpath'));
addpath(fullfile(directory, 'matlab'));
defaults = sd_experiment_config();
ordering = defaults.defaultOrdering;
config = sd_experiment_config(mode, ordering);
paths = setup_sd_experiment_paths();
manifest = prepare_matrix_cache(mode, "", ~defaults.useCachedData, family, ordering);
settings = manifest(string(manifest.Family) == family, :);
if family == "random-dimension"
    settings = settings(ismember(settings.Dimension, config.dimensionValues), :);
end
settings = sortrows(settings, {'Dimension', 'TargetCondition'});
summary = table();
records = table();
names = {'Dimension', 'ConditionNumber', 'SampleIndex', 'CDFLOPs', ...
    'SDFLOPs', 'CDIGTs', 'SDIGTs', 'CDSwaps', 'SDSwaps', ...
    'SDInheritanceIGTs', 'SDReductionIGTs', 'SDInheritanceSwaps', ...
    'SDReductionSwaps'};

for settingIndex = 1:height(settings)
    cached = load(char(string(settings.File(settingIndex))), 'matrices');
    n = settings.Dimension(settingIndex);
    count = size(cached.matrices, 3);
    values = zeros(count, numel(names));
    for sampleIndex = 1:count
        covariance = cached.matrices(:, :, sampleIndex);
        child = covariance(1:n-1, 1:n-1);
        parent = lambda4_reduce_covariance_counted(covariance);
        [cd, cdFlops, cdEvents] = lambda4_reduce_covariance_counted(child);
        [sd, sdFlops, sdEvents] = sd_lambda4_reduce_child_counted(covariance, parent);
        cdResidual = norm(cd.Qz - cd.Z * child * cd.Z', 'fro');
        sdResidual = norm(sd.Qz - sd.Z * child * sd.Z', 'fro');
        assert(max(cdResidual, sdResidual) / max(1, norm(child, 'fro')) <= 1e-4, ...
            'A counted reduction failed the covariance consistency check.');
        values(sampleIndex, :) = [n, settings.TargetCondition(settingIndex), ...
            sampleIndex, cdFlops, sdFlops, cdEvents.integerGaussianTransforms, ...
            sdEvents.totalIGT, cdEvents.conditionalSwaps, sdEvents.totalSwaps, ...
            sdEvents.inheritanceIGT, sdEvents.reductionIGT, ...
            sdEvents.inheritanceSwaps, sdEvents.reductionSwaps];
    end
    records = [records; array2table(values, 'VariableNames', names)]; %#ok<AGROW>
    row = array2table([n, settings.TargetCondition(settingIndex), count, ...
        mean(values(:, 4:end), 1)], 'VariableNames', ...
        [names(1:2), {'NumMatrices'}, names(4:end)]);
    summary = [summary; row]; %#ok<AGROW>
    fprintf('%s: n=%d, %d matrices counted.\n', family, n, count);
end

output = fullfile(paths.numericalOutput, mode, ordering, 'operation_counts');
if ~isfolder(output), mkdir(output); end
writetable(summary, fullfile(output, family + "_summary.csv"));
writetable(records, fullfile(output, family + "_per_matrix.csv"));
disp(summary);
end
