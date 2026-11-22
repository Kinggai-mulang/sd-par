function results = run_gb_chain_comparison(mode)
%RUN_GB_CHAIN_COMPARISON Generate inputs and compare complete GB ILS sequences.
arguments
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])} = "pilot"
end

directory = fileparts(mfilename('fullpath'));
addpath(fullfile(directory, 'matlab'));
defaults = sd_experiment_config();
ordering = defaults.defaultOrdering;
config = sd_experiment_config(mode, ordering);
setup_sd_experiment_paths();
prepare_matrix_cache(mode, "", ~defaults.useCachedData, "geometry-based", ordering);
results = run_deletion_chain_experiment(mode, config.timingRepetitions, ...
    "", ordering, config.numRandomMatrices, "geometry-based");
disp(results);
end
