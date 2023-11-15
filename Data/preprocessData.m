clear, clc, close all

cd('/Volumes/bc6/PROJECTS/CausalGlobs/analysis')

cleandata = analysis('dirname',{'210422_135417_Pollux', '211012_124119_Pollux'}, 'heatmapfig', true, 'makefigs', false);
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/210422--211012_Pollux_data.mat','cleandata')

cleandata = analysis('dirname',{'210517_090630_Castor', '211108_090705_Castor'}, 'heatmapfig', true, 'makefigs', false);
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/210517--211108_Castor_data.mat','cleandata')

cleandata = analysis('dirname',{'210428_125912_Buster', '210609_124628_Buster'}, 'heatmapfig', true, 'makefigs', false);
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/210428--210609_Buster_data.mat','cleandata')

cleandata = analysis('dirname',{'220322_091149_Morty', '220823_081207_Morty'}, 'heatmapfig', true, 'makefigs', false);
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/220322--220823_Morty_data.mat','cleandata')


%% Testing getting more variables...

% cleandata = analysis('dirname',{'210422_135417_Pollux'}, 'heatmapfig', true, 'choicefixfig', true, 'locationfig', true, 'makefigs', false);
