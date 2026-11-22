function [Qa, metadata] = generate_controlled_spd( ...
        dimension, conditionNumber, seed, ambiguityOrdering)
%GENERATE_CONTROLLED_SPD Random SPD matrix with prescribed 2-norm condition.
%   The eigenvalues form a geometric progression and are normalized so
%   trace(QA) equals DIMENSION. Ambiguities are then ordered by increasing
%   marginal variance, so the last ambiguity is removed first.

arguments
    dimension (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(dimension, 2)}
    conditionNumber (1, 1) double {mustBeFinite, mustBeGreaterThanOrEqual(conditionNumber, 1)}
    seed (1, 1) double {mustBeInteger, mustBeNonnegative}
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

stream = RandStream('mt19937ar', 'Seed', seed);
[orthogonalMatrix, triangularMatrix] = qr(randn(stream, dimension));
columnSigns = sign(diag(triangularMatrix));
columnSigns(columnSigns == 0) = 1;
orthogonalMatrix = orthogonalMatrix .* columnSigns';

if dimension == 2
    exponents = [0 1];
else
    exponents = (0:dimension-1) / (dimension-1);
end
eigenvalues = conditionNumber .^ exponents;
eigenvalues = dimension * eigenvalues / sum(eigenvalues);

Qa = orthogonalMatrix * diag(eigenvalues) * orthogonalMatrix';
Qa = (Qa + Qa') / 2;
if ambiguityOrdering == "variance-ordered"
    [~, order] = sort(diag(Qa), 'ascend');
else
    order = randperm(stream, dimension);
end
Qa = Qa(order, order);

metadata.generator = 'controlled-spd';
metadata.dimension = dimension;
metadata.targetCondition = conditionNumber;
metadata.actualCondition = cond(Qa, 2);
metadata.trace = trace(Qa);
metadata.seed = seed;
metadata.permutation = order;
metadata.ordering = char(ambiguityOrdering);
metadata.removal = 'remove last after ordering';

validate_covariance_matrix( ...
    Qa, dimension, conditionNumber, ambiguityOrdering);

end
