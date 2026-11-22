function timing = benchmark_cd_sd_reduction(Qparent, repetitions)
%BENCHMARK_CD_SD_REDUCTION Time uninstrumented CD and SD reductions.
%   FLOP counters are deliberately not called here. CD and SD execution
%   order alternates between repetitions to limit order-dependent bias.

arguments
    Qparent (:, :) double {mustBeReal, mustBeFinite}
    repetitions (1, 1) double {mustBeInteger, mustBePositive} = 20
end

n = size(Qparent, 1);
if size(Qparent, 2) ~= n || n < 2
    error('benchmark_cd_sd_reduction:InvalidDimensions', ...
        'Qparent must be a square matrix of dimension at least two.');
end

parent = lambda4_reduce_covariance(Qparent);
Qchild = Qparent(1:n-1, 1:n-1);
checkedChild = sd_lambda4_reduce_child(Qparent, parent, n);
normalizer = max(1, norm(Qchild, 'fro'));
residual = norm(checkedChild.Qz ...
    - checkedChild.Z * Qchild * checkedChild.Z', 'fro') / normalizer;
if residual > 1e-4
    error('benchmark_cd_sd_reduction:CovarianceResidual', ...
        'SD covariance residual %.3g exceeds the tolerance.', residual);
end

% Exercise both complete paths before either timer starts. Two passes remove
% one-time JIT, toolbox-loading, and dimension-dependent allocation costs.
for warmupIndex = 1:2
    lambda4_reduce_covariance(Qchild);
    sd_lambda4_reduce_child(Qparent, parent, n);
end

% Measure several operations inside each timed batch. This makes timer
% overhead negligible at low dimensions while alternating CD/SD order to
% avoid a systematic first- or second-run advantage.
numBatches = min(4, repetitions);
batchSizes = floor(repetitions / numBatches) * ones(numBatches, 1);
batchSizes(1:mod(repetitions, numBatches)) = ...
    batchSizes(1:mod(repetitions, numBatches)) + 1;
cdBatchSeconds = zeros(numBatches, 1);
sdBatchSeconds = zeros(numBatches, 1);
for batchIndex = 1:numBatches
    if mod(batchIndex, 2) == 1
        startTime = tic;
        for repetition = 1:batchSizes(batchIndex)
            lambda4_reduce_covariance(Qchild);
        end
        cdBatchSeconds(batchIndex) = toc(startTime);

        startTime = tic;
        for repetition = 1:batchSizes(batchIndex)
            sd_lambda4_reduce_child(Qparent, parent, n);
        end
        sdBatchSeconds(batchIndex) = toc(startTime);
    else
        startTime = tic;
        for repetition = 1:batchSizes(batchIndex)
            sd_lambda4_reduce_child(Qparent, parent, n);
        end
        sdBatchSeconds(batchIndex) = toc(startTime);

        startTime = tic;
        for repetition = 1:batchSizes(batchIndex)
            lambda4_reduce_covariance(Qchild);
        end
        cdBatchSeconds(batchIndex) = toc(startTime);
    end
end
cdSingle = cdBatchSeconds ./ batchSizes;
sdSingle = sdBatchSeconds ./ batchSizes;

timing.repetitions = repetitions;
timing.cdSeconds = cdSingle;
timing.sdSeconds = sdSingle;
timing.cdMeanSeconds = sum(cdBatchSeconds) / repetitions;
timing.sdMeanSeconds = sum(sdBatchSeconds) / repetitions;
timing.cdMedianSeconds = median(cdSingle);
timing.sdMedianSeconds = median(sdSingle);

end
