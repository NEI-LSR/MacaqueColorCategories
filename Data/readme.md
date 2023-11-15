## Raw Data
Raw data available in the Causal Globs repo (`CausalGlobs/data`).

## Preprocessed data
Data extracted to here by running `preprocessData.m` which calls `CausalGlobs/Analyses/analysis.m` 

TODO Make it so that it incorporates cue/choice locations, and trial timing

## Conversion to csv
Converting to csv to make it accessible to folks without MATLAB (and also smaller filesize).

Using: `MacaqueColorCategories/Data/mat2csv.m`

## Combining data

With `combineData.m` for csv files.
With `combineData_mat.m` for mat files.