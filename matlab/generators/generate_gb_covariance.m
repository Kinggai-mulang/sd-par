function [Qa, metadata] = generate_gb_covariance( ...
        settings, numSatellites, seed, ambiguityOrdering)
%GENERATE_GB_COVARIANCE Triple-frequency short-baseline GPS GB covariance.
%   Extends the single-frequency design in LAMBDA 4.0 RUN_example_3.
%   Azimuths and elevations are independent idealized sky draws.

arguments
    settings (1, 1) struct
    numSatellites (1, 1) double {mustBeInteger, mustBePositive}
    seed (1, 1) double {mustBeInteger, mustBeNonnegative}
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])}
end

if settings.numReceivers ~= 2 || settings.numEpochs ~= 1 || ...
        settings.numFrequencies ~= numel(settings.frequenciesHz)
    error('generate_gb_covariance:InvalidSettings', ...
        'The GB generator requires two receivers, one epoch and one frequency per entry.');
end
stream = RandStream('mt19937ar', 'Seed', seed);
elevation = (settings.minElevationDegrees + ...
    (settings.maxElevationDegrees-settings.minElevationDegrees) ...
    * rand(stream, numSatellites, 1)) * pi/180;
azimuth = (2*rand(stream, numSatellites, 1)-1) * pi;
lineOfSight = [cos(elevation).*sin(azimuth), ...
    cos(elevation).*cos(azimuth), sin(elevation)];
numDD = numSatellites-1;
DD = [-ones(numDD, 1), eye(numDD)];
G = DD*lineOfSight;
if rank(G) < 3
    error('generate_gb_covariance:SingularGeometry', ...
        'The sampled sky geometry does not determine the three baseline components.');
end

J = settings.numFrequencies;
wavelength = 299792458 ./ settings.frequenciesHz(:);
H = [repmat(G, J, 1), zeros(J*numDD, J*numDD); ...
    repmat(G, J, 1), kron(diag(wavelength), eye(numDD))];
Qdd = DD*DD';
Qcode = 2*settings.codeStdMeters^2*Qdd;
Qphase = 2*settings.phaseStdMeters^2*Qdd;
Qy = blkdiag(kron(eye(J), Qcode), kron(eye(J), Qphase));
normal = H'*(Qy\H);
if rcond(normal) < 1e-14
    error('generate_gb_covariance:SingularGeometry', ...
        'The GB normal matrix is numerically singular.');
end
Qx = normal\eye(size(normal));
Qa = Qx(4:end, 4:end);
Qa = (Qa+Qa')/2;
if ambiguityOrdering == "variance-ordered"
    [~, permutation] = sort(diag(Qa), 'ascend');
else
    orderStream = RandStream('mt19937ar', 'Seed', seed+500000000);
    permutation = randperm(orderStream, size(Qa, 1));
end
Qa = Qa(permutation, permutation);
validate_covariance_matrix(Qa, J*numDD, [], ambiguityOrdering);

metadata.actualCondition = cond(Qa, 2);
metadata.permutation = permutation;
metadata.azimuthRadians = azimuth;
metadata.elevationRadians = elevation;
metadata.seed = seed;
end
