function comparison = compare_lambda4_child_reduction(Qparent, removeIndex, repetitions)
%COMPARE_LAMBDA4_CHILD_REDUCTION Prepare one conventional-versus-SD case.
%   COMPARISON = COMPARE_LAMBDA4_CHILD_REDUCTION(QPARENT, REMOVEINDEX,
%   REPETITIONS) times conventional child-set reduction and sequential
%   decorrelation with the same LAMBDA 4.0 implementation. Parent reduction
%   is performed once and excluded from both child-set timings.

arguments
    Qparent (:, :) double {mustBeReal, mustBeFinite}
    removeIndex (1, 1) double {mustBeInteger, mustBePositive}
    repetitions (1, 1) double {mustBeInteger, mustBePositive} = 100
end

parent = lambda4_reduce_covariance(Qparent);
inherited = sd_inherit_child_transform(Qparent, parent, removeIndex);

conventionalElapsed = zeros(repetitions, 1);
sequentialElapsed = zeros(repetitions, 1);
for repetition = 1:repetitions
    startTime = tic;
    conventional = lambda4_reduce_covariance(inherited.Qsubset);
    conventionalElapsed(repetition) = toc(startTime);

    startTime = tic;
    sequential = sd_lambda4_reduce_child(Qparent, parent, removeIndex);
    sequentialElapsed(repetition) = toc(startTime);
end

conventionalModel = conventional.Z * inherited.Qsubset * conventional.Z';
sequentialModel = sequential.Z * inherited.Qsubset * sequential.Z';
normalizer = max(1, norm(inherited.Qsubset, 'fro'));

comparison.dimension = size(inherited.Qsubset, 1);
comparison.removeIndex = removeIndex;
comparison.repetitions = repetitions;
comparison.conventionalTimeSeconds = mean(conventionalElapsed);
comparison.sequentialTimeSeconds = mean(sequentialElapsed);
comparison.speedup = comparison.conventionalTimeSeconds / ...
    comparison.sequentialTimeSeconds;
comparison.conventionalResidual = ...
    norm(conventional.Qz - conventionalModel, 'fro') / normalizer;
comparison.sequentialResidual = ...
    norm(sequential.Qz - sequentialModel, 'fro') / normalizer;
comparison.conventional = conventional;
comparison.sequential = sequential;

end
