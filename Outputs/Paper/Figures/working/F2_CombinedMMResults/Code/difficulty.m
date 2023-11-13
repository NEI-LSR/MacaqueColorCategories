clear, clc, close all

% csv or mat?
csv = 0;

%%

filename = 'combinedData';
repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',];
addpath(genpath([repoHomeDir,filesep,'Analyses']))

if csv
    warning('Code for using csv is not complete')

    loadedData = readtable([repoHomeDir,filesep,'Data',filesep,'combinedData.csv']);

    for i = 1:size(loadedData,1)
        data.trialdata.cues{i,1} = table2array(loadedData(i,2));
        data.trialdata.choices{i,1} = table2array(loadedData(i,4:7));
        data.trialdata.chosen{i,1} = table2array(loadedData(i,3));
    end
else
    load([repoHomeDir,filesep,'Data',filesep,'combinedData.mat']);
end

%%

difficulty_psychometric(cleandata,filename);



