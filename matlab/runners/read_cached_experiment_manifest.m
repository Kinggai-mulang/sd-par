function manifest = read_cached_experiment_manifest( ...
        cacheRoot, mode, ambiguityOrdering, useCachedData)
%READ_CACHED_EXPERIMENT_MANIFEST Read an existing cache manifest only.

arguments
    cacheRoot (1, 1) string
    mode (1, 1) string {mustBeMember(mode, ["pilot", "full"])}
    ambiguityOrdering (1, 1) string {mustBeMember(ambiguityOrdering, ...
        ["variance-ordered", "random-order"])}
    useCachedData (1, 1) logical
end

if ~useCachedData
    error('read_cached_experiment_manifest:CacheReuseRequired', ...
        'This experiment reads an existing matrix cache and never rebuilds it.');
end

if ambiguityOrdering == "random-order"
    cacheRoot = fullfile(cacheRoot, 'random_order');
end
manifestFile = fullfile(cacheRoot, "manifest_" + mode + ".tsv");
if ~exist(manifestFile, 'file')
    error('read_cached_experiment_manifest:MissingManifest', ...
        'No existing %s cache manifest was found.', mode);
end
manifest = readtable(manifestFile, 'FileType', 'text', ...
    'Delimiter', '\t', 'VariableNamingRule', 'preserve');

end
