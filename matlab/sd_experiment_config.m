function config = sd_experiment_config(mode, ambiguityOrdering)
%SD_EXPERIMENT_CONFIG User-editable settings for every SD experiment.

arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])} = "variance-ordered"
end

% ------------------------- Editable settings --------------------------
config.experimentSwitches = logical([0 1 0 0 0 0 0]);
config.experimentMode = "pilot";
config.useCachedData = true;
config.useCachedResults = true; % Experiment 7 checkpoints only.
config.defaultOrdering = "random-order";
config.lambdaDirectory = ""; % Set an external LAMBDA 4.0 directory if needed.
config.stageDimensionValues = 10:1:50;  % Experiment 1 (includes search).
config.dimensionValues = 10:2:80;       % Experiment 2 (reduction only).
config.conditionDimension = 40;
config.conditionExponents = 2:2:20; % Tests kappa = 2.^exponents.
config.numRandomMatricesPilot = 3;
config.numRandomMatricesFull = 1000;
config.timingRepetitionsPilot = 2;
config.timingRepetitionsFull = 20;
config.gbSearchBatchRepetitionsPilot = 3;
config.gbSearchBatchRepetitionsFull = 10;
config.baseSeed = 20260921;
config.ambiguityOrdering = ambiguityOrdering;

% Experiment 4: independent idealized GPS sky geometries (one float
% ambiguity sample per covariance), following the GB design of LAMBDA's
% RUN_example_3 and extending it to L1/L2/L5.
config.gb.numFrequencies = 3;
config.gb.numReceivers = 2;
config.gb.numEpochs = 1;
config.gb.numSatellites = 5:3:32;
config.gb.frequenciesHz = [1575.42e6; 1227.60e6; 1176.45e6];
config.gb.codeStdMeters = 0.300;
config.gb.phaseStdMeters = 0.003;
config.gb.minElevationDegrees = 15;
config.gb.maxElevationDegrees = 75;
% ---------------------------------------------------------------------

config.mode = mode;
if mode == "pilot"
    config.numRandomMatrices = config.numRandomMatricesPilot;
    config.timingRepetitions = config.timingRepetitionsPilot;
    config.gbSearchBatchRepetitions = config.gbSearchBatchRepetitionsPilot;
else
    config.numRandomMatrices = config.numRandomMatricesFull;
    config.timingRepetitions = config.timingRepetitionsFull;
    config.gbSearchBatchRepetitions = config.gbSearchBatchRepetitionsFull;
end

validateattributes(config.stageDimensionValues, {'double'}, ...
    {'row','integer','>=',2,'increasing'});
validateattributes(config.dimensionValues, {'double'}, ...
    {'row','integer','>=',2,'increasing'});
validateattributes(config.conditionDimension, {'double'}, ...
    {'scalar','integer','>=',2});
validateattributes(config.conditionExponents, {'double'}, ...
    {'row','integer','positive','increasing'});
validateattributes(config.numRandomMatrices, {'double'}, ...
    {'scalar','integer','positive'});
validateattributes(config.timingRepetitions, {'double'}, ...
    {'scalar','integer','positive'});
validateattributes(config.gbSearchBatchRepetitions, {'double'}, ...
    {'scalar','integer','positive'});
end
