function cleandata = processBaeData()

%% Load Bae data

% clear, clc, close all

% TODO How did we go from the csvs to the .mat file?

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
filepath = [repoHomeDir,'Data',filesep,'secondaryData',filesep,'Bae_2015',filesep,'processed_data',filesep];

filename = 'bae_human_data.mat';

load([filepath,filename])

%%

% target color = cue index
% cue angle = angle of target color on each trial
% choice angle = angle of choice on each trial

colresp = [cue_angle, choice_angle];

all_cues = length(unique(target_color)); % number of cue colors in original data
interval = 360/all_cues; % degrees between each cue
cues_degrees = target_color*interval;
cues = target_color;

cues_degrees(cues_degrees==360) = 0;

error = zeros(size(cue_angle));
for i = 1:size(colresp,1)% for trial
    if abs(colresp(i,2) - colresp(i,1)) < 180
        error(i,1) = (colresp(i,2) - colresp(i,1));
    elseif abs(colresp(i,2) - colresp(i,1)) > 180 && colresp(i,2) < colresp(i,1)
        error(i,1) = ((360 - colresp(i,1) + colresp(i,2)));
    else
        error(i,1) = (abs(colresp(i,2) - colresp(i,1)) - 360);
    end
end

chosen_degrees = cues_degrees + error;
chosen_degrees(chosen_degrees > 360 - interval) = chosen_degrees(chosen_degrees > 360 - interval) - 360;
chosen_degrees(chosen_degrees < 0) = chosen_degrees(chosen_degrees < 0) + 360;

chosen_ind = round(chosen_degrees/interval) + 1;
% chosen_ind = round(chosen_degrees/interval);

choices = repmat({1:180},length(cues),1);

nTrials = size(colresp,1);

cleandata.trialdata.cues = num2cell(cues);
cleandata.trialdata.chosen = num2cell(chosen_ind);
cleandata.trialdata.choices = choices;
cleandata.trialdata.dirname = repmat({'Bae_Humans'},nTrials,1);
cleandata.trialdata.paradigm = repmat({'Bae_Humans'},nTrials,1);
cleandata.trialdata.allchoices = repmat({1:1:180},nTrials,1);

end
