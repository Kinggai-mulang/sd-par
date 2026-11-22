function result = sd_lambda4_reduce_child(Qparent, parentReduction, removeIndex)
%SD_LAMBDA4_REDUCE_CHILD Sequential decorrelation using LAMBDA 4.0.
%   The parent Qz and inverse Z are updated incrementally during inheritance;
%   the final reduction uses the unmodified LAMBDA 4.0 core functions.

inherited = sd_inherit_child_transform( ...
    Qparent, parentReduction, removeIndex);
finalReduction = lambda4_reduce_inherited_covariance( ...
    inherited.Qinitial, inherited.inverseBsubset);

result.Qsubset = inherited.Qsubset;
result.Bsubset = inherited.Bsubset;
result.Qinitial = inherited.Qinitial;
result.Z = finalReduction.Z;
result.Qz = finalReduction.Qz;
result.L = finalReduction.L;
result.d = finalReduction.d;
result.inverseZ = finalReduction.inverseZ;
result.permutation = inherited.permutation;

end
