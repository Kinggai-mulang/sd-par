function result = lambda4_reduce_covariance(Qa)
%LAMBDA4_REDUCE_COVARIANCE Conventional LAMBDA 4.0 reduction baseline.
%   RESULT = LAMBDA4_REDUCE_COVARIANCE(QA) reduces the ambiguity
%   variance-covariance matrix QA with decorrelateVC from LAMBDA 4.0.
%   RESULT.Z uses the manuscript convention QZ = Z*QA*Z'.

arguments
    Qa (:, :) double {mustBeReal, mustBeFinite}
end

if size(Qa, 1) ~= size(Qa, 2)
    error('lambda4_reduce_covariance:NonSquareInput', ...
        'Qa must be a square variance-covariance matrix.');
end
if exist('decorrelateVC', 'file') ~= 2
    error('lambda4_reduce_covariance:MissingLambda4', ...
        ['LAMBDA 4.0 is not on the MATLAB path. Add the repository-level ', ...
         'external LAMBDA_toolbox directory before calling this function.']);
end

% decorrelateVC only assembles its explicit Z matrix when an ambiguity
% vector is supplied. A zero vector requests that output without changing
% the covariance reduction.
[Qz, L, d, inverseZ, toolboxZ] = decorrelateVC(Qa, zeros(size(Qa, 1), 1));

% LAMBDA 4.0 defines z = toolboxZ' * a. The manuscript instead writes the
% same left-acting transform as Qz = Z * Qa * Z'.
result.Z = round(toolboxZ');
result.Qz = Qz;
result.L = L;
result.d = d;
result.inverseZ = round(inverseZ);

end
