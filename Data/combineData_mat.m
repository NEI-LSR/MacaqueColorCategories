function cleandata = combineData_mat(dirname,rn,bootstrap)

if ~exist('rn','var')
    rn = 2;
    warning('Using default random number')
end
rng(rn); % fix random number generator for reproducibility

if ~exist('bootstrap','var')
    boostrap = false;
end

%% Load data

filename{1} = '210422--211012_Pollux_data';
filename{2} = '210517--211108_Castor_data.mat';
filename{3} = '210428--210609_Buster_data';
filename{4} = '220322--220823_Morty_data.mat';

for dataset = 1:length(filename)
    tempdata{dataset} = load([dirname,filesep,filename{dataset}]);
end

%% Strip out aborts

for i = 1:length(tempdata)
    abortIndex = isnan(cell2mat(tempdata{1,i}.trialdata.chosen_idx));
    tempdata{1,i}.trialdata.allchoices  = tempdata{1,i}.trialdata.allchoices(~abortIndex);
    tempdata{1,i}.trialdata.dirname     = tempdata{1,i}.trialdata.dirname(~abortIndex);
    tempdata{1,i}.trialdata.paradigm    = tempdata{1,i}.trialdata.paradigm(~abortIndex);
    tempdata{1,i}.trialdata.choices     = tempdata{1,i}.trialdata.choices(~abortIndex);
    tempdata{1,i}.trialdata.chosen_idx  = tempdata{1,i}.trialdata.chosen_idx(~abortIndex);
    tempdata{1,i}.trialdata.cues        = tempdata{1,i}.trialdata.cues(~abortIndex);
    tempdata{1,i}.trialdata.stimCols_raw= tempdata{1,i}.trialdata.stimCols_raw(~abortIndex);
    tempdata{1,i}.trialdata.chosen      = tempdata{1,i}.trialdata.chosen(~abortIndex);
    tempdata{1,i}.trialdata.stimCols    = tempdata{1,i}.trialdata.stimCols(~abortIndex);
end

%%

for dataset = 1:length(filename)
    nTrials(dataset) = length(tempdata{1,dataset}.trialdata.cues);
end

for dataset = 1:length(filename)
    if ~boostrap
        idx(dataset,:) = randsample(nTrials(dataset),min(nTrials));
    elseif bootstrap
        idx(dataset,:) = randi(nTrials(dataset),min(nTrials),1);
    end
end

%% Subset

cleandata.trialdata = struct();

cleandata.trialdata.choices =...
    [   tempdata{1,1}.trialdata.choices(idx(1,:));...
        tempdata{1,2}.trialdata.choices(idx(2,:));...
        tempdata{1,3}.trialdata.choices(idx(3,:));...
        tempdata{1,4}.trialdata.choices(idx(4,:))];

cleandata.trialdata.chosen_idx =...
    [   tempdata{1,1}.trialdata.chosen_idx(idx(1,:));...
        tempdata{1,2}.trialdata.chosen_idx(idx(2,:));...
        tempdata{1,3}.trialdata.chosen_idx(idx(3,:));...
        tempdata{1,4}.trialdata.chosen_idx(idx(4,:))];

cleandata.trialdata.cues =...
    [   tempdata{1,1}.trialdata.cues(idx(1,:));...
        tempdata{1,2}.trialdata.cues(idx(2,:));...
        tempdata{1,3}.trialdata.cues(idx(3,:));...
        tempdata{1,4}.trialdata.cues(idx(4,:))];

cleandata.trialdata.chosen =...
    [   tempdata{1,1}.trialdata.chosen(idx(1,:));...
        tempdata{1,2}.trialdata.chosen(idx(2,:));...
        tempdata{1,3}.trialdata.chosen(idx(3,:));...
        tempdata{1,4}.trialdata.chosen(idx(4,:))];

cleandata.trialdata.stimCols =...
    [   tempdata{1,1}.trialdata.stimCols(idx(1,:));...
        tempdata{1,2}.trialdata.stimCols(idx(2,:));...
        tempdata{1,3}.trialdata.stimCols(idx(3,:));...
        tempdata{1,4}.trialdata.stimCols(idx(4,:))];

cleandata.trialdata.paradigm =...
    [   tempdata{1,1}.trialdata.paradigm(idx(1,:));...
        tempdata{1,2}.trialdata.paradigm(idx(2,:));...
        tempdata{1,3}.trialdata.paradigm(idx(3,:));...
        tempdata{1,4}.trialdata.paradigm(idx(4,:))];

cleandata.trialdata.dirname = {'combined'};

%% Whole set

% cleandata.trialdata = struct();
% 
% cleandata.trialdata.choices =...
%     [   tempdata{1,1}.trialdata.choices;...
%         tempdata{1,2}.trialdata.choices;...
%         tempdata{1,3}.trialdata.choices;...
%         tempdata{1,4}.trialdata.choices];
% 
% cleandata.trialdata.chosen_idx =...
%     [   tempdata{1,1}.trialdata.chosen_idx;...
%         tempdata{1,2}.trialdata.chosen_idx;...
%         tempdata{1,3}.trialdata.chosen_idx;...
%         tempdata{1,4}.trialdata.chosen_idx];
% 
% cleandata.trialdata.cues =...
%     [   tempdata{1,1}.trialdata.cues;...
%         tempdata{1,2}.trialdata.cues;...
%         tempdata{1,3}.trialdata.cues;...
%         tempdata{1,4}.trialdata.cues];
% 
% cleandata.trialdata.chosen =...
%     [   tempdata{1,1}.trialdata.chosen;...
%         tempdata{1,2}.trialdata.chosen;...
%         tempdata{1,3}.trialdata.chosen;...
%         tempdata{1,4}.trialdata.chosen];
% 
% cleandata.trialdata.stimCols =...
%     [   tempdata{1,1}.trialdata.stimCols,...
%         tempdata{1,2}.trialdata.stimCols,...
%         tempdata{1,3}.trialdata.stimCols,...
%         tempdata{1,4}.trialdata.stimCols];
% 
% cleandata.trialdata.paradigm =...
%     [   tempdata{1,1}.trialdata.paradigm,...
%         tempdata{1,2}.trialdata.paradigm,...
%         tempdata{1,3}.trialdata.paradigm,...
%         tempdata{1,4}.trialdata.paradigm];
% 
% cleandata.trialdata.dirname = {'combined'};

%% Subset of each animal individually

% i = 4;
% 
% cleandata.trialdata = struct();
% 
% cleandata.trialdata.choices =...
%     [   tempdata{1,i}.trialdata.choices(idx(i,:))];
% 
% cleandata.trialdata.chosen_idx =...
%     [   tempdata{1,i}.trialdata.chosen_idx(idx(i,:))];
% 
% cleandata.trialdata.cues =...
%     [   tempdata{1,i}.trialdata.cues(idx(i,:))];
% 
% cleandata.trialdata.chosen =...
%     [   tempdata{1,i}.trialdata.chosen(idx(i,:))];
% 
% cleandata.trialdata.stimCols =...
%     [   tempdata{1,i}.trialdata.stimCols(idx(i,:))];
% 
% cleandata.trialdata.paradigm =...
%     [   tempdata{1,i}.trialdata.paradigm(idx(i,:))];
% 
% cleandata.trialdata.dirname = {'combined'};

end