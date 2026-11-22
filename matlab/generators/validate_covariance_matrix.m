function diagnostics = validate_covariance_matrix( ...
        Qa, expectedDimension, expectedCondition, ambiguityOrdering)
%VALIDATE_COVARIANCE_MATRIX Fail fast on invalid cached experiment input.

arguments
    Qa (:, :) double {mustBeReal, mustBeFinite}
    expectedDimension (1, 1) double {mustBeInteger, mustBePositive}
    expectedCondition double = []
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

if ~isequal(size(Qa), [expectedDimension expectedDimension])
    error('validate_covariance_matrix:WrongDimension', ...
        'Expected a %d-by-%d covariance matrix.', ...
        expectedDimension, expectedDimension);
end

normalizer = max(1, norm(Qa, 'fro'));
symmetryResidual = norm(Qa - Qa', 'fro') / normalizer;
if symmetryResidual > 1e-12
    error('validate_covariance_matrix:NotSymmetric', ...
        'Relative symmetry residual %.3g exceeds the tolerance.', ...
        symmetryResidual);
end

[~, positiveDefiniteFlag] = chol(Qa);
if positiveDefiniteFlag ~= 0
    error('validate_covariance_matrix:NotPositiveDefinite', ...
        'The covariance matrix is not numerically positive definite.');
end

if exist('decomposeLtDL', 'file') == 2
    try
        decomposeLtDL(Qa);
    catch cause
        error('validate_covariance_matrix:LambdaDecompositionFailure', ...
            'LAMBDA 4.0 rejected the covariance: %s', cause.message);
    end
end

if ambiguityOrdering == "variance-ordered" && ...
        any(diff(diag(Qa)) < -1e-12 * max(1, max(diag(Qa))))
    error('validate_covariance_matrix:WrongOrdering', ...
        'Marginal variances must be nondecreasing.');
end

actualCondition = cond(Qa, 2);
if ~isempty(expectedCondition) && ~isnan(expectedCondition)
    relativeConditionError = abs(actualCondition - expectedCondition) / ...
        max(1, expectedCondition);
    if relativeConditionError > 1e-8
        error('validate_covariance_matrix:WrongCondition', ...
            'Relative condition-number error %.3g exceeds the tolerance.', ...
            relativeConditionError);
    end
else
    relativeConditionError = NaN;
end

diagnostics.symmetryResidual = symmetryResidual;
diagnostics.actualCondition = actualCondition;
diagnostics.relativeConditionError = relativeConditionError;
diagnostics.minimumCholeskyDiagonal = min(diag(chol(Qa)));

end
