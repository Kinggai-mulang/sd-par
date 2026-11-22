function [result, flops, operations] = ...
        sd_lambda4_reduce_child_counted(Qparent, parentReduction)
%SD_LAMBDA4_REDUCE_CHILD_COUNTED Count one last-ambiguity SD reduction.
%   Inherited covariance and inverse transformations are updated directly;
%   no dense Bsub*Qsubset*Bsub' product is formed.

arguments
    Qparent (:, :) double {mustBeReal, mustBeFinite}
    parentReduction (1, 1) struct
end

n = size(Qparent, 1);
requiredFields = {'Z', 'Qz', 'inverseZ'};
if n < 2 || size(Qparent, 2) ~= n || ...
        ~all(isfield(parentReduction, requiredFields))
    error('sd_lambda4_reduce_child_counted:InvalidParent', ...
        'A square covariance and complete parent reduction are required.');
end

Qsubset = Qparent(1:n-1, 1:n-1);
B = round(parentReduction.Z);
inverseB = round(parentReduction.inverseZ);
transformedCovariance = parentReduction.Qz;
flops = 0;
operations.inheritanceIGT = 0;
operations.inheritanceSwaps = 0;

while true
    column = B(:, n);
    nonzeroRows = find(column ~= 0);
    if isscalar(nonzeroRows) && abs(column(nonzeroRows)) == 1
        break
    end
    if isempty(nonzeroRows)
        error('sd_lambda4_reduce_child_counted:SingularTransform', ...
            'The inherited last column vanished.');
    end

    [~, minimumPosition] = min(abs(column(nonzeroRows)));
    pivotRow = nonzeroRows(minimumPosition);
    pivotValue = B(pivotRow, n);
    for row = reshape(nonzeroRows, 1, [])
        if row == pivotRow
            continue
        end
        multiplier = round(B(row, n) / pivotValue);
        B(row, :) = B(row, :) - multiplier * B(pivotRow, :);
        transformedCovariance(row, :) = transformedCovariance(row, :) ...
            - multiplier * transformedCovariance(pivotRow, :);
        transformedCovariance(:, row) = transformedCovariance(:, row) ...
            - multiplier * transformedCovariance(:, pivotRow);
        inverseB(:, pivotRow) = inverseB(:, pivotRow) ...
            + multiplier * inverseB(:, row);
        flops = flops + 1 + 8 * n;
        if multiplier ~= 0
            operations.inheritanceIGT = operations.inheritanceIGT + 1;
        end
    end
end

unitRow = find(B(:, n) ~= 0, 1);
if unitRow ~= n
    B([unitRow n], :) = B([n unitRow], :);
    transformedCovariance([unitRow n], :) = ...
        transformedCovariance([n unitRow], :);
    transformedCovariance(:, [unitRow n]) = ...
        transformedCovariance(:, [n unitRow]);
    inverseB(:, [unitRow n]) = inverseB(:, [n unitRow]);
    operations.inheritanceSwaps = operations.inheritanceSwaps + 1;
end
if B(n, n) == -1
    B(n, :) = -B(n, :);
    transformedCovariance(n, :) = -transformedCovariance(n, :);
    transformedCovariance(:, n) = -transformedCovariance(:, n);
    inverseB(:, n) = -inverseB(:, n);
    flops = flops + 4 * n;
end
for row = 1:n-1
    multiplier = B(row, n);
    if multiplier == 0
        continue
    end
    B(row, :) = B(row, :) - multiplier * B(n, :);
    transformedCovariance(row, :) = transformedCovariance(row, :) ...
        - multiplier * transformedCovariance(n, :);
    transformedCovariance(:, row) = transformedCovariance(:, row) ...
        - multiplier * transformedCovariance(:, n);
    inverseB(:, n) = inverseB(:, n) + multiplier * inverseB(:, row);
    flops = flops + 8 * n;
    operations.inheritanceIGT = operations.inheritanceIGT + 1;
end

B = round(B);
inverseB = round(inverseB);
Bsubset = B(1:n-1, 1:n-1);
inverseBsubset = inverseB(1:n-1, 1:n-1);
Qinitial = transformedCovariance(1:n-1, 1:n-1);

[finalReduction, finalFlops, reductionOperations] = ...
    lambda4_reduce_covariance_counted(Qinitial, inverseBsubset);
flops = flops + finalFlops;
operations.reductionIGT = reductionOperations.integerGaussianTransforms;
operations.reductionSwaps = reductionOperations.conditionalSwaps;
operations.totalIGT = operations.inheritanceIGT + operations.reductionIGT;
operations.totalSwaps = operations.inheritanceSwaps + operations.reductionSwaps;

result.Qsubset = Qsubset;
result.Bsubset = Bsubset;
result.Qinitial = Qinitial;
result.Z = finalReduction.Z;
result.Qz = finalReduction.Qz;
result.L = finalReduction.L;
result.d = finalReduction.d;
result.inverseZ = finalReduction.inverseZ;

end
