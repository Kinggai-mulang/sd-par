function [result, flops, operations] = ...
        lambda4_reduce_covariance_counted(Qa, inverseInitial)
%LAMBDA4_REDUCE_COVARIANCE_COUNTED LAMBDA 4.0 reduction with FLOP count.
%   Algorithm source: LAMBDA 4.0 by Lotfi Massarweh, June 2024.
%   Source copyright: Geoscience & Remote Sensing department, TU Delft.
%   See THIRD_PARTY_NOTICES.md for attribution and dependency information.
%   This reference implementation follows decomposeLtDL, transformZ, and
%   computeIGT_row from LAMBDA 4.0. Addition, subtraction, multiplication,
%   and division each count as one FLOP. Comparisons, indexing, ROUND, and
%   data movement do not count. The explicit inverse used to recover Z is
%   assigned the standard LU-plus-solves cost 8*n^3/3.
%
%   This function is for operation counting only. Runtime experiments must
%   call the uninstrumented lambda4_reduce_covariance implementation.

arguments
    Qa (:, :) double {mustBeReal, mustBeFinite}
    inverseInitial double = []
end

n = size(Qa, 1);
if n < 1 || size(Qa, 2) ~= n
    error('lambda4_reduce_covariance_counted:InvalidCovariance', ...
        'Qa must be a nonempty square matrix.');
end

flops = 0;
operations.integerGaussianTransforms = 0;
operations.conditionalSwaps = 0;
working = Qa;
for k = n:-1:1
    m = k - 1;
    if m == 0
        continue
    end
    working(k, 1:m) = working(k, 1:m) / working(k, k);
    flops = flops + m;
    working(1:m, 1:m) = working(1:m, 1:m) ...
        - working(k, 1:m)' * working(k, k) * working(k, 1:m);
    flops = flops + m + m^2 + m^2;
end

L = tril(working, -1) + eye(n);
d = diag(working)';
if any(d < 1e-12)
    error('lambda4_reduce_covariance_counted:NotPositiveDefinite', ...
        'LAMBDA 4.0 positive-definiteness threshold was not satisfied.');
end

if isempty(inverseInitial)
    inverseZ = eye(n);
elseif isequal(size(inverseInitial), [n n])
    inverseZ = inverseInitial;
else
    error('lambda4_reduce_covariance_counted:DimensionMismatch', ...
        'inverseInitial must be empty or have the same dimension as Qa.');
end
k = n - 1;
while k > 0
    kp1 = k + 1;
    correlation = L(kp1, k);
    mu = round(correlation);
    if mu ~= 0
        correlation = correlation - mu;
        flops = flops + 1;
    end

    delta = d(k) + correlation^2 * d(kp1);
    flops = flops + 3;
    if delta < d(kp1)
        if mu ~= 0
            rowCount = n - kp1 + 1;
            L(kp1:n, k) = L(kp1:n, k) - mu * L(kp1:n, kp1);
            inverseZ(:, kp1) = inverseZ(:, kp1) + mu * inverseZ(:, k);
            flops = flops + 2 * rowCount + 2 * n;
            operations.integerGaussianTransforms = ...
                operations.integerGaussianTransforms + 1;

            for row = k+2:n
                mu = round(L(row, k));
                if mu ~= 0
                    rowCount = n - row + 1;
                    L(row:n, k) = L(row:n, k) - mu * L(row:n, row);
                    inverseZ(:, row) = inverseZ(:, row) + mu * inverseZ(:, k);
                    flops = flops + 2 * rowCount + 2 * n;
                    operations.integerGaussianTransforms = ...
                        operations.integerGaussianTransforms + 1;
                end
            end
        end

        lambda = L(kp1, k) * d(kp1) / delta;
        eta = d(k) / delta;
        flops = flops + 3;

        columnCount = k - 1;
        swapMatrix = [-L(kp1, k), 1; eta, lambda];
        L([k kp1], 1:columnCount) = ...
            swapMatrix * L([k kp1], 1:columnCount);
        flops = flops + 6 * columnCount;

        L(kp1, k) = lambda;
        d(k) = eta * d(kp1);
        d(kp1) = delta;
        flops = flops + 1;

        L(k+2:n, [k kp1]) = L(k+2:n, [kp1 k]);
        inverseZ(:, [k kp1]) = inverseZ(:, [kp1 k]);
        if k < n - 1
            k = k + 1;
        end
        operations.conditionalSwaps = operations.conditionalSwaps + 1;
    else
        k = k - 1;
    end
end

for row = 2:n
    muValues = round(L(row, 1:row-1));
    nonzeroColumns = find(muValues);
    for column = nonzeroColumns
        rowCount = n - row + 1;
        L(row:n, column) = L(row:n, column) ...
            - muValues(column) * L(row:n, row);
        inverseZ(:, row) = inverseZ(:, row) ...
            + muValues(column) * inverseZ(:, column);
        flops = flops + 2 * rowCount + 2 * n;
        operations.integerGaussianTransforms = ...
            operations.integerGaussianTransforms + 1;
    end
end

scaledL = d' .* L;
Qz = L' * scaledL;
flops = flops + n^2 + n^3 + n^2 * (n - 1);

toolboxZ = round(inv(inverseZ'));
flops = flops + (8 / 3) * n^3;

result.Z = round(toolboxZ');
result.Qz = Qz;
result.L = L;
result.d = d;
result.inverseZ = round(inverseZ);

end
