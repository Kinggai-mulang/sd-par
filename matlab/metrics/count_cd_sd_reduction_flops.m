function counts = count_cd_sd_reduction_flops(Qparent)
%COUNT_CD_SD_REDUCTION_FLOPS Count CD and SD for one child problem.

arguments
    Qparent (:, :) double {mustBeReal, mustBeFinite}
end

n = size(Qparent, 1);
if size(Qparent, 2) ~= n || n < 2
    error('count_cd_sd_reduction_flops:InvalidDimensions', ...
        'Qparent must be a square matrix of dimension at least two.');
end

[parent, ~] = lambda4_reduce_covariance_counted(Qparent);
Qchild = Qparent(1:n-1, 1:n-1);
[cdChild, counts.cd] = lambda4_reduce_covariance_counted(Qchild);
[sdChild, counts.sd] = ...
    sd_lambda4_reduce_child_counted(Qparent, parent);
assertEquivalent(Qchild, cdChild);
assertEquivalent(Qchild, sdChild);

counts.speedup = counts.cd / counts.sd;

end

function assertEquivalent(Qa, reduction)
normalizer = max(1, norm(Qa, 'fro'));
residual = norm(reduction.Qz - reduction.Z * Qa * reduction.Z', 'fro') ...
    / normalizer;
if residual > 1e-4
    error('count_cd_sd_reduction_flops:ResidualTooLarge', ...
        'Normalized covariance residual %.3g exceeds the tolerance.', residual);
end
end
