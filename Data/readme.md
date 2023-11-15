## Raw Data
Raw data available in the Causal Globs repo (`CausalGlobs/data`).

## Preprocessed data
Data extracted to here by running `CausalGlobs/Analyses/analysis.m` as:

```
cleandata = analysis('dirname',{'210422_135417_Pollux', '211012_124119_Pollux'}, 'heatmapfig', true, 'makefigs', false);

cleandata = analysis('dirname',{'210517_090630_Castor', '211108_090705_Castor'}, 'heatmapfig', true, 'makefigs', false)

cleandata = analysis('dirname',{'210428_125912_Buster', '210609_124628_Buster'}, 'heatmapfig', true, 'makefigs', false)

cleandata = analysis('dirname',{'220322_091149_Morty', '220823_081207_Morty'}, 'heatmapfig', true, 'makefigs', false)
```

TODO Make it so that it incorporates cue/choice locations, and trial timing

## Conversion to csv
Converting to csv to make it accessible to folks without MATLAB (and also smaller filesize).

Using: `MacaqueColorCategories/Data/mat2csv.m`

## Combining data

With `combineData.m` for csv files.
With `combineData_mat.m` for mat files.