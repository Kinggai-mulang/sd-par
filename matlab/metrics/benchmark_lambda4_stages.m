function timing = benchmark_lambda4_stages(Qa, aHat, repetitions)
%BENCHMARK_LAMBDA4_STAGES Time LAMBDA 4.0 reduction and ILS search alone.
%   This is the Experiment 1 primitive. It calls no SD implementation.

arguments
    Qa (:, :) double {mustBeReal, mustBeFinite}
    aHat (:, 1) double {mustBeReal, mustBeFinite}
    repetitions (1, 1) double {mustBeInteger, mustBePositive} = 20
end

n = size(Qa, 1);
if size(Qa, 2) ~= n || numel(aHat) ~= n
    error('benchmark_lambda4_stages:DimensionMismatch', ...
        'Qa and aHat must describe the same ambiguity problem.');
end
if exist('estimatorILS', 'file') ~= 2
    error('benchmark_lambda4_stages:MissingLambda4', ...
        'Run setup_sd_experiment_paths before benchmarking.');
end

reduction = lambda4_reduce_covariance(Qa);
zHat = reduction.Z * aHat;

% Exercise both complete code paths before starting either timer. Repeating
% the warm-up avoids charging MATLAB JIT compilation, lazy toolbox loading,
% and the first allocation of dimension-dependent work arrays to the first
% (usually lowest-dimensional) measurement.
for warmupIndex = 1:2
    lambda4_reduce_covariance(Qa);
    estimatorILS(zHat, reduction.L, reduction.d, 2);
end

% Time the short reductions as one batch. Starting and stopping a timer for
% every call otherwise biases the lowest dimensions, where timer overhead
% is comparable with the operation itself.
startTime = tic;
for repetition = 1:repetitions
    lambda4_reduce_covariance(Qa);
end
reductionBatchSeconds = toc(startTime);

% Choose an adaptive batch size from an untimed-result pilot batch. Very
% short searches are repeated until the formal timing window is about 5 ms;
% intrinsically hard searches remain single calls. The cap provides a firm
% upper bound on the extra work.
pilotSearchRepetitions = 3;
startTime = tic;
for repetition = 1:pilotSearchRepetitions
    estimatorILS(zHat, reduction.L, reduction.d, 2);
end
pilotSearchSeconds = toc(startTime);
searchRepetitions = min(1000, max(1, ceil( ...
    0.005 * pilotSearchRepetitions / pilotSearchSeconds)));

startTime = tic;
for repetition = 1:searchRepetitions
    estimatorILS(zHat, reduction.L, reduction.d, 2);
end
searchBatchSeconds = toc(startTime);
searchTimes = repmat( ...
    searchBatchSeconds / searchRepetitions, searchRepetitions, 1);
reductionTimes = repmat(reductionBatchSeconds / repetitions, repetitions, 1);

timing.repetitions = repetitions;
timing.reductionSeconds = reductionTimes;
timing.searchSeconds = searchTimes;
timing.reductionMeanSeconds = mean(reductionTimes);
timing.searchMeanSeconds = mean(searchTimes);
timing.reductionMedianSeconds = median(reductionTimes);
timing.searchMedianSeconds = median(searchTimes);

end
