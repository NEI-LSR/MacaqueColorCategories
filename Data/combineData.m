function combinedData = combineData(DataDir)

rng(2); % fix random number generator for reproducibility

%% Load data

filename{1} = '210422--211012_Pollux_data';
filename{2} = '210517--211108_Castor_data';
filename{3} = '220322--220823_Morty_data';
filename{4} = '210428--210609_Buster_data';

for participant = 1:length(filename)
    tempdata{participant} = readtable([DataDir,filesep,filename{participant},'.csv']);
end

%% Strip out aborts

for participant = 1:length(tempdata)

    abortIndex = isnan(tempdata{1,participant}.chosen) | ... 
        any(isnan([...
        tempdata{1,participant}.choices_1,...
        tempdata{1,participant}.choices_2,...
        tempdata{1,participant}.choices_3,...
        tempdata{1,participant}.choices_4]),2);

    tempdata{1,participant} = tempdata{1,participant}(~abortIndex,:);

end

%% Subsetting

% Which dataset has fewest trials?
for dataset = 1:length(filename)
    nTrials(dataset) = length(tempdata{1,dataset}.cues);
end

% Generate random indexes to subsample all datasets down to match the one
% with fewest trials.
for dataset = 1:length(filename)
    idx(dataset,:) = randsample(nTrials(dataset),min(nTrials));
end

%% Subset, matching contribution from each participant

combinedData = [...
    tempdata{1,1}(idx(1,:),:);...
    tempdata{1,2}(idx(2,:),:);...
    tempdata{1,3}(idx(3,:),:);...
    tempdata{1,4}(idx(4,:),:)];

%% Whole set

% combinedData = [...
%     tempdata{1,1};...
%     tempdata{1,2};...
%     tempdata{1,3};...
%     tempdata{1,4}];

%% Subset of each animal individually

% participant = 2;
% combinedData = tempdata{1,participant}(idx(participant,:),:);

%%

writetable(combinedData,[DataDir,filesep,'combinedData.csv']);

end