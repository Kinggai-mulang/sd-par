function [Qa, metadata] = generate_random_ldlt(dimension, seed, ambiguityOrdering)
%GENERATE_RANDOM_LDLT Random covariance for the dimension experiments.
%   D has independent diagonal entries in [0.2, 1.2]. Strictly lower
%   entries of the unit lower-triangular L are independent in [-1, 1].
%   Qa=L'*D*L follows the original dimension-simulation model.

arguments
    dimension (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(dimension, 2)}
    seed (1, 1) double {mustBeInteger, mustBeNonnegative}
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

stream = RandStream('mt19937ar', 'Seed', seed);
L = eye(dimension) + tril(2 * rand(stream, dimension) - 1, -1);
d = 0.2 + rand(stream, dimension, 1);
Qa = L' * (d .* L);
Qa = (Qa + Qa') / 2;
if ambiguityOrdering == "variance-ordered"
    [~, order] = sort(diag(Qa), 'ascend');
else
    order = randperm(stream, dimension);
end
Qa = Qa(order, order);

metadata.generator = 'random-ldlt';
metadata.dimension = dimension;
metadata.actualCondition = cond(Qa, 2);
metadata.seed = seed;
metadata.diagonalRange = [0.2 1.2];
metadata.lowerTriangularRange = [-1 1];
metadata.permutation = order;
metadata.ordering = char(ambiguityOrdering);
metadata.removal = 'remove last after ordering';

validate_covariance_matrix(Qa, dimension, [], ambiguityOrdering);

end
