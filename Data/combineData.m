function cleandata = combineData(dirname,rn,bootstrap,eccentricity)

% To combine data from individual files into one
% matching number of trials so that each monkey is equally represented
%
% Standard usage:
% > cleandata = combineData('.')
% > save("combinedData.mat","cleandata")

if ~exist('rn','var') || isempty(rn)
    rn = 5; % random number seed picked by hand
    warning('Using default random number')
end
rng(rn); % fix random number generator for reproducibility

if ~exist('bootstrap','var') || isempty(bootstrap)
    bootstrap = false;
end

%% Load data

filename{1} = '210422--211012_Pollux_data.mat';
filename{2} = '210517--211108_Castor_data.mat';
filename{3} = '210428--210609_Buster_data.mat';
filename{4} = '220322--220823_Morty_data.mat';

for dataset = 1:length(filename)
    tempdata{dataset} = load([dirname,filesep,filename{dataset}]);
end

%% Subset by eccentricity

if exist('eccentricity','var')
    for dataset = 1:length(filename)
        for trial = 1:length(tempdata{dataset}.cleandata.trialdata.dirname)
            cue_x = tempdata{1, dataset}.cleandata.trialdata.cue_loc{trial, 1}(1);
            cue_y = tempdata{1, dataset}.cleandata.trialdata.cue_loc{trial, 1}(2);
            tempdata{1,dataset}.cleandata.trialdata.cue_eccentricity(trial,1) = sqrt((cue_x-(1920/2))^2 + (cue_y-(1080/2))^2);
        end
        if strcmp(eccentricity,'upper') %remove certain eccentricities
            eccIndex = tempdata{dataset}.cleandata.trialdata.cue_eccentricity(:,1) >= mean(tempdata{dataset}.cleandata.trialdata.cue_eccentricity(:,1));
        else
            eccIndex = tempdata{dataset}.cleandata.trialdata.cue_eccentricity(:,1) < mean(tempdata{dataset}.cleandata.trialdata.cue_eccentricity(:,1));
        end
        tempdata{1,dataset}.cleandata.trialdata.dirname     = tempdata{1,dataset}.cleandata.trialdata.dirname(eccIndex);
        tempdata{1,dataset}.cleandata.trialdata.paradigm    = tempdata{1,dataset}.cleandata.trialdata.paradigm(eccIndex);
        tempdata{1,dataset}.cleandata.trialdata.choices     = tempdata{1,dataset}.cleandata.trialdata.choices(eccIndex);
        tempdata{1,dataset}.cleandata.trialdata.cues        = tempdata{1,dataset}.cleandata.trialdata.cues(eccIndex);
        tempdata{1,dataset}.cleandata.trialdata.chosen      = tempdata{1,dataset}.cleandata.trialdata.chosen(eccIndex);
    end
end

%% Strip out aborts

for i = 1:length(tempdata)
    nTrialsAll(i) = size(tempdata{1,i}.cleandata.trialdata.chosen,1);
    nAborts(i) = sum(isnan(cell2mat(tempdata{1,i}.cleandata.trialdata.chosen)));
    abortRate(i) = nAborts(i)/nTrialsAll(i);
end

for i = 1:length(tempdata)
    abortIndex = isnan(cell2mat(tempdata{1,i}.cleandata.trialdata.chosen));
    tempdata{1,i}.cleandata.trialdata.dirname     = tempdata{1,i}.cleandata.trialdata.dirname(~abortIndex);
    tempdata{1,i}.cleandata.trialdata.paradigm    = tempdata{1,i}.cleandata.trialdata.paradigm(~abortIndex);
    tempdata{1,i}.cleandata.trialdata.choices     = tempdata{1,i}.cleandata.trialdata.choices(~abortIndex);
    tempdata{1,i}.cleandata.trialdata.cues        = tempdata{1,i}.cleandata.trialdata.cues(~abortIndex);
    tempdata{1,i}.cleandata.trialdata.chosen      = tempdata{1,i}.cleandata.trialdata.chosen(~abortIndex);
end

%%

for dataset = 1:length(filename)
    nTrials(dataset) = length(tempdata{1,dataset}.cleandata.trialdata.cues);
end

for dataset = 1:length(filename)
    if ~bootstrap
        idx(dataset,:) = randsample(nTrials(dataset),min(nTrials));
    elseif bootstrap
        idx(dataset,:) = randi(nTrials(dataset),min(nTrials),1);
    end
end

%% Subset

cleandata.trialdata = struct();

cleandata.trialdata.choices =...
    [   tempdata{1,1}.cleandata.trialdata.choices(idx(1,:));...
        tempdata{1,2}.cleandata.trialdata.choices(idx(2,:));...
        tempdata{1,3}.cleandata.trialdata.choices(idx(3,:));...
        tempdata{1,4}.cleandata.trialdata.choices(idx(4,:))];

cleandata.trialdata.cues =...
    [   tempdata{1,1}.cleandata.trialdata.cues(idx(1,:));...
        tempdata{1,2}.cleandata.trialdata.cues(idx(2,:));...
        tempdata{1,3}.cleandata.trialdata.cues(idx(3,:));...
        tempdata{1,4}.cleandata.trialdata.cues(idx(4,:))];

cleandata.trialdata.chosen =...
    [   tempdata{1,1}.cleandata.trialdata.chosen(idx(1,:));...
        tempdata{1,2}.cleandata.trialdata.chosen(idx(2,:));...
        tempdata{1,3}.cleandata.trialdata.chosen(idx(3,:));...
        tempdata{1,4}.cleandata.trialdata.chosen(idx(4,:))];

cleandata.trialdata.paradigm =...
    [   tempdata{1,1}.cleandata.trialdata.paradigm(idx(1,:)),...
        tempdata{1,2}.cleandata.trialdata.paradigm(idx(2,:)),...
        tempdata{1,3}.cleandata.trialdata.paradigm(idx(3,:)),...
        tempdata{1,4}.cleandata.trialdata.paradigm(idx(4,:))];

cleandata.trialdata.dirname = cell(size(cleandata.trialdata.cues'));
cleandata.trialdata.dirname(:) = {'combined'};


%% Whole set

% cleandata.trialdata = struct();
% 
% cleandata.cleandata.trialdata.choices =...
%     [   tempdata{1,1}.cleandata.trialdata.choices;...
%         tempdata{1,2}.cleandata.trialdata.choices;...
%         tempdata{1,3}.cleandata.trialdata.choices;...
%         tempdata{1,4}.cleandata.trialdata.choices];
% 
% cleandata.cleandata.trialdata.chosen_idx =...
%     [   tempdata{1,1}.cleandata.trialdata.chosen_idx;...
%         tempdata{1,2}.cleandata.trialdata.chosen_idx;...
%         tempdata{1,3}.cleandata.trialdata.chosen_idx;...
%         tempdata{1,4}.cleandata.trialdata.chosen_idx];
% 
% cleandata.cleandata.trialdata.cues =...
%     [   tempdata{1,1}.cleandata.trialdata.cues;...
%         tempdata{1,2}.cleandata.trialdata.cues;...
%         tempdata{1,3}.cleandata.trialdata.cues;...
%         tempdata{1,4}.cleandata.trialdata.cues];
% 
% cleandata.cleandata.trialdata.chosen =...
%     [   tempdata{1,1}.cleandata.trialdata.chosen;...
%         tempdata{1,2}.cleandata.trialdata.chosen;...
%         tempdata{1,3}.cleandata.trialdata.chosen;...
%         tempdata{1,4}.cleandata.trialdata.chosen];
% 
% cleandata.cleandata.trialdata.stimCols =...
%     [   tempdata{1,1}.cleandata.trialdata.stimCols,...
%         tempdata{1,2}.cleandata.trialdata.stimCols,...
%         tempdata{1,3}.cleandata.trialdata.stimCols,...
%         tempdata{1,4}.cleandata.trialdata.stimCols];
% 
% cleandata.cleandata.trialdata.paradigm =...
%     [   tempdata{1,1}.cleandata.trialdata.paradigm,...
%         tempdata{1,2}.cleandata.trialdata.paradigm,...
%         tempdata{1,3}.cleandata.trialdata.paradigm,...
%         tempdata{1,4}.cleandata.trialdata.paradigm];
% 
% cleandata.cleandata.trialdata.dirname = {'combined'};

%% Subset of each animal individually

% i = 4;
% 
% cleandata.trialdata = struct();
% 
% cleandata.cleandata.trialdata.choices =...
%     [   tempdata{1,i}.cleandata.trialdata.choices(idx(i,:))];
% 
% cleandata.cleandata.trialdata.chosen_idx =...
%     [   tempdata{1,i}.cleandata.trialdata.chosen_idx(idx(i,:))];
% 
% cleandata.cleandata.trialdata.cues =...
%     [   tempdata{1,i}.cleandata.trialdata.cues(idx(i,:))];
% 
% cleandata.cleandata.trialdata.chosen =...
%     [   tempdata{1,i}.cleandata.trialdata.chosen(idx(i,:))];
% 
% cleandata.cleandata.trialdata.stimCols =...
%     [   tempdata{1,i}.cleandata.trialdata.stimCols(idx(i,:))];
% 
% cleandata.cleandata.trialdata.paradigm =...
%     [   tempdata{1,i}.cleandata.trialdata.paradigm(idx(i,:))];
% 
% cleandata.cleandata.trialdata.dirname = {'combined'};

end