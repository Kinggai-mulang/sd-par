function inherited = sd_inherit_child_transform(Qparent, parentReduction, removeIndex)
%SD_INHERIT_CHILD_TRANSFORM Inherit a child transform after one exclusion.
%   PARENTREDUCTION is the result of lambda4_reduce_covariance. Integer row
%   operations are applied simultaneously to Z, Qz, and inverse(Z), so the
%   inherited covariance is obtained without forming Bsub*Qsubset*Bsub'.

arguments
    Qparent (:, :) double {mustBeReal, mustBeFinite}
    parentReduction (1, 1) struct
    removeIndex (1, 1) double {mustBeInteger, mustBePositive}
end

n = size(Qparent, 1);
requiredFields = {'Z', 'Qz', 'inverseZ'};
if n < 2 || size(Qparent, 2) ~= n
    error('sd_inherit_child_transform:InvalidCovariance', ...
        'Qparent must be a square matrix of dimension at least two.');
end
if ~all(isfield(parentReduction, requiredFields))
    error('sd_inherit_child_transform:IncompleteParentReduction', ...
        'parentReduction must contain Z, Qz, and inverseZ.');
end
if removeIndex > n
    error('sd_inherit_child_transform:InvalidIndex', ...
        'removeIndex must not exceed the parent dimension.');
end

Zparent = parentReduction.Z;
if ~isequal(size(Zparent), [n n]) || ...
        ~isequal(size(parentReduction.Qz), [n n]) || ...
        ~isequal(size(parentReduction.inverseZ), [n n])
    error('sd_inherit_child_transform:DimensionMismatch', ...
        'The parent reduction must have the same dimension as Qparent.');
end
integerTolerance = 1e-9;
if any(abs(Zparent(:) - round(Zparent(:))) > integerTolerance)
    error('sd_inherit_child_transform:NonIntegerTransform', ...
        'The parent Z transformation must be integer valued.');
end

order = [1:removeIndex-1, removeIndex+1:n, removeIndex];
permutation = eye(n);
permutation = permutation(order, :);
Qpermuted = Qparent(order, order);
Qsubset = Qpermuted(1:n-1, 1:n-1);

% B = Zparent*permutation'. Therefore inv(B) = permutation*inv(Zparent).
B = round(Zparent(:, order));
inverseB = round(parentReduction.inverseZ(order, :));
transformedCovariance = parentReduction.Qz;
lastColumn = n;

while true
    column = B(:, lastColumn);
    nonzeroRows = find(column ~= 0);
    if isscalar(nonzeroRows) && abs(column(nonzeroRows)) == 1
        break
    end
    if isempty(nonzeroRows)
        error('sd_inherit_child_transform:SingularTransform', ...
            'The last column vanished during integer row reduction.');
    end

    [~, minimumPosition] = min(abs(column(nonzeroRows)));
    pivotRow = nonzeroRows(minimumPosition);
    pivotValue = B(pivotRow, lastColumn);
    for row = reshape(nonzeroRows, 1, [])
        if row == pivotRow
            continue
        end
        multiplier = round(B(row, lastColumn) / pivotValue);
        B(row, :) = B(row, :) - multiplier * B(pivotRow, :);
        transformedCovariance(row, :) = transformedCovariance(row, :) ...
            - multiplier * transformedCovariance(pivotRow, :);
        transformedCovariance(:, row) = transformedCovariance(:, row) ...
            - multiplier * transformedCovariance(:, pivotRow);
        inverseB(:, pivotRow) = inverseB(:, pivotRow) ...
            + multiplier * inverseB(:, row);
    end
end

unitRow = find(B(:, lastColumn) ~= 0, 1);
if unitRow ~= n
    B([unitRow n], :) = B([n unitRow], :);
    transformedCovariance([unitRow n], :) = ...
        transformedCovariance([n unitRow], :);
    transformedCovariance(:, [unitRow n]) = ...
        transformedCovariance(:, [n unitRow]);
    inverseB(:, [unitRow n]) = inverseB(:, [n unitRow]);
end
if B(n, lastColumn) == -1
    B(n, :) = -B(n, :);
    transformedCovariance(n, :) = -transformedCovariance(n, :);
    transformedCovariance(:, n) = -transformedCovariance(:, n);
    inverseB(:, n) = -inverseB(:, n);
end
for row = 1:n-1
    multiplier = B(row, lastColumn);
    if multiplier == 0
        continue
    end
    B(row, :) = B(row, :) - multiplier * B(n, :);
    transformedCovariance(row, :) = transformedCovariance(row, :) ...
        - multiplier * transformedCovariance(n, :);
    transformedCovariance(:, row) = transformedCovariance(:, row) ...
        - multiplier * transformedCovariance(:, n);
    inverseB(:, n) = inverseB(:, n) + multiplier * inverseB(:, row);
end

B = round(B);
inverseB = round(inverseB);
if ~isequal(B(:, lastColumn), [zeros(n-1, 1); 1])
    error('sd_inherit_child_transform:ReductionFailure', ...
        'The inherited transformation did not reach the required block form.');
end

inherited.Qsubset = Qsubset;
inherited.Bsubset = B(1:n-1, 1:n-1);
inherited.inverseBsubset = inverseB(1:n-1, 1:n-1);
inherited.Qinitial = transformedCovariance(1:n-1, 1:n-1);
inherited.permutation = permutation;
inherited.blockTransform = B;

end
