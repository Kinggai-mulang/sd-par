function paths = setup_sd_experiment_paths()
%SETUP_SD_EXPERIMENT_PATHS Add only the directories required by the study.

matlabDirectory = fileparts(mfilename('fullpath'));
experimentDirectory = fileparts(matlabDirectory);
config = sd_experiment_config();
lambdaDirectory = char(config.lambdaDirectory);
if isempty(lambdaDirectory)
    lambdaDirectory = getenv('SD_LAMBDA_DIRECTORY');
end
if isempty(lambdaDirectory)
    error('setup_sd_experiment_paths:MissingToolbox', ...
        'Set config.lambdaDirectory or the SD_LAMBDA_DIRECTORY environment variable.');
end
if ~exist(fullfile(lambdaDirectory, 'LAMBDA_toolbox', 'decorrelateVC.m'), 'file')
    error('setup_sd_experiment_paths:MissingToolbox', ...
        'Set config.lambdaDirectory to your LAMBDA 4.0 folder.');
end

paths.matlab = matlabDirectory;
paths.experiment = experimentDirectory;
paths.lambda = lambdaDirectory;
paths.cache = fullfile(experimentDirectory, 'data', 'generated');
paths.output = fullfile(experimentDirectory, 'output');
paths.numericalOutput = fullfile(paths.output, 'numerical');
paths.figureOutput = fullfile(paths.output, 'figures');

addpath(fullfile(matlabDirectory, 'algorithms'));
addpath(fullfile(matlabDirectory, 'generators'));
addpath(fullfile(matlabDirectory, 'metrics'));
addpath(fullfile(matlabDirectory, 'plotting'));
addpath(fullfile(matlabDirectory, 'runners'));
addpath(fullfile(lambdaDirectory, 'LAMBDA_toolbox'));
if isfolder(fullfile(lambdaDirectory, 'LAMBDA_examples'))
    addpath(fullfile(lambdaDirectory, 'LAMBDA_examples'));
end

end
