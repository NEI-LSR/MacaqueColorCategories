clear, clc, close all

cd('/Volumes/bc6/PROJECTS/CausalGlobs/analysis')
% cd('D:\20230203\CausalGlobs\analysis')

cleandata = analysis('dirname',{'210422_135417_Pollux', '211012_124119_Pollux'},...
    'colCats_dataExtract', true);
cleandata.trialdata = rmfield(cleandata.trialdata,{'choice_onset','chosen_idx','holdtoselect_ms','absolute_RT','trial_length'}); % removing unrequired fields to reduce file size
cleandata = rmfield(cleandata,'timestampdata');
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/210422--211012_Pollux_data.mat','cleandata')
% save('C:\Users\cege-user\Documents\MacaqueColorCategories\Data\210422--211012_Pollux_data.mat','cleandata')

cleandata = analysis('dirname',{'210517_090630_Castor', '211108_090705_Castor'},... 
    'colCats_dataExtract', true);
cleandata.trialdata = rmfield(cleandata.trialdata,{'choice_onset','chosen_idx','holdtoselect_ms','absolute_RT','trial_length'}); % removing unrequired fields to reduce file size
cleandata = rmfield(cleandata,'timestampdata');
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/210517--211108_Castor_data.mat','cleandata')
% save('C:\Users\cege-user\Documents\MacaqueColorCategories\Data\210517--211108_Castor_data.mat','cleandata')

cleandata = analysis('dirname',{'210428_125912_Buster', '210609_124628_Buster'},...
    'colCats_dataExtract', true);
cleandata.trialdata = rmfield(cleandata.trialdata,{'choice_onset','chosen_idx','holdtoselect_ms','absolute_RT','trial_length'}); % removing unrequired fields to reduce file size
cleandata = rmfield(cleandata,'timestampdata');
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/210428--210609_Buster_data.mat','cleandata')
% save('C:\Users\cege-user\Documents\MacaqueColorCategories\Data\210428--210609_Buster_data.mat','cleandata')

cleandata = analysis('dirname',{'220322_091149_Morty', '220823_081207_Morty'},...
    'colCats_dataExtract', true);
cleandata.trialdata = rmfield(cleandata.trialdata,{'choice_onset','chosen_idx','holdtoselect_ms','absolute_RT','trial_length'}); % removing unrequired fields to reduce file size
cleandata = rmfield(cleandata,'timestampdata');
save('/Volumes/bc6/PROJECTS/MacaqueColorCategories/Data/220322--220823_Morty_data.mat','cleandata')
% save('C:\Users\cege-user\Documents\MacaqueColorCategories\Data\220322--220823_Morty_data.mat','cleandata')


