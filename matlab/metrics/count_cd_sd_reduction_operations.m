function operations = count_cd_sd_reduction_operations(Qparent)
%COUNT_CD_SD_REDUCTION_OPERATIONS Count integer reduction events for one child.
%   Counts nonzero integer Gauss row operations and actual swaps. For SD,
%   inheritance operations and final LAMBDA reduction operations are kept
%   separate and then added for the total comparison.

arguments
    Qparent (:, :) double {mustBeReal, mustBeFinite}
end

n = size(Qparent, 1);
if size(Qparent, 2) ~= n || n < 2
    error('count_cd_sd_reduction_operations:InvalidDimensions', ...
        'Qparent must be a square matrix of dimension at least two.');
end

[parentReduction, ~] = lambda4_reduce_covariance_counted(Qparent);
Qchild = Qparent(1:n-1, 1:n-1);
[cdReduction, ~, cdEvents] = ...
    lambda4_reduce_covariance_counted(Qchild);
[sdReduction, ~, sdEvents] = ...
    sd_lambda4_reduce_child_counted(Qparent, parentReduction);

assertEquivalent(Qchild, cdReduction);
assertEquivalent(Qchild, sdReduction);

operations.cdReductionIGT = cdEvents.integerGaussianTransforms;
operations.sdInheritanceIGT = sdEvents.inheritanceIGT;
operations.sdReductionIGT = sdEvents.reductionIGT;
operations.sdTotalIGT = sdEvents.totalIGT;
operations.cdReductionSwaps = cdEvents.conditionalSwaps;
operations.sdInheritanceSwaps = sdEvents.inheritanceSwaps;
operations.sdReductionSwaps = sdEvents.reductionSwaps;
operations.sdTotalSwaps = sdEvents.totalSwaps;

end

function assertEquivalent(Qa, reduction)
normalizer = max(1, norm(Qa, 'fro'));
residual = norm(reduction.Qz - reduction.Z * Qa * reduction.Z', 'fro') ...
    / normalizer;
if residual > 1e-4
    error('count_cd_sd_reduction_operations:ResidualTooLarge', ...
        'Normalized covariance residual %.3g exceeds the tolerance.', residual);
end
end
