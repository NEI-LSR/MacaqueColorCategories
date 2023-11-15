clear, clc, close all

% csv or mat?
csv = 0;

%% sets paths

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

difficulty_psychometric(data,filename);

%% to get the number of cues per animal that have been subsampled to match Buster

cues = cell2mat(data.trialdata.cues);

for i = 1:64
    t(i,1) = sum(cues(1:24526) == i);
end

for i = 1:64
    t(i,2) = sum(cues(24527:49052) == i);
end

for i = 1:64
    t(i,3) = sum(cues(49053:73578) == i);
end

for i = 1:64
    t(i,4) = sum(cues(73579:end) == i);
end

disp(sum(t(:)) == size(cues,1)) % sanity check



