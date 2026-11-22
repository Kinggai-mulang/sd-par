function result = lambda4_reduce_inherited_covariance(Qinitial, inverseInitial)
%LAMBDA4_REDUCE_INHERITED_COVARIANCE Finish reduction from an inherited Z.
%   Uses the unmodified LAMBDA 4.0 decomposition and transform functions.
%   INVERSEINITIAL is inv(Bsubset) for Qinitial=Bsubset*Qsubset*Bsubset'.

arguments
    Qinitial (:, :) double {mustBeReal, mustBeFinite}
    inverseInitial (:, :) double {mustBeReal, mustBeFinite}
end

n = size(Qinitial, 1);
if size(Qinitial, 2) ~= n || ~isequal(size(inverseInitial), [n n])
    error('lambda4_reduce_inherited_covariance:DimensionMismatch', ...
        'Qinitial and inverseInitial must be equally sized square matrices.');
end
if exist('decomposeLtDL', 'file') ~= 2 || exist('transformZ', 'file') ~= 2
    error('lambda4_reduce_inherited_covariance:MissingLambda4', ...
        'Run setup_sd_experiment_paths before reducing a child covariance.');
end

[L, d] = decomposeLtDL(Qinitial);
[L, d, inverseZ] = transformZ(L, d, inverseInitial);
Qz = L' * (d' .* L);
Z = round(inv(inverseZ));

result.Z = Z;
result.Qz = Qz;
result.L = L;
result.d = d;
result.inverseZ = round(inverseZ);

end
