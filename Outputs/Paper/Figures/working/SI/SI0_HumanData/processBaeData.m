function cleandata = processBaeData()

% Load Bae data
% filepath = 'C:\Users\selwynhm\Documents\MATLAB\nAFC\analysis\';
%filepath = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\protocol\categories\Reference_Data\Bae_2015\processed_data\';
% filepath = 'C:\Users\garsidedj\Documents\CausalGlobs\protocol\categories\Reference_Data\Bae_2015\processed_data\';
filepath = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\protocol\categories\Reference_Data\Bae_2015\processed_data\';

filename = 'bae_human_data.mat';

load([filepath,filename])

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

chosen = cues_degrees + error;
chosen(chosen > 360) = chosen(chosen > 360) - 360;
chosen(chosen < 0) = chosen(chosen < 0) + 360;

chosen = chosen/2; % chosen = chosen/2;
chosen = round(chosen);

chosen(chosen == 0) = 1;


choices = repmat({1:180},length(cues),1);

% cues = cues/4;
% cues = round(cues);
% cues(cues == 0) = 1;

nTrials = size(colresp,1);

cleandata.trialdata.cues = num2cell(cues);
cleandata.trialdata.chosen = num2cell(chosen);
cleandata.trialdata.choices = choices;
cleandata.trialdata.dirname = repmat({'Bae_Humans'},nTrials,1);
cleandata.trialdata.paradigm = repmat({'Bae_Humans'},nTrials,1);
cleandata.trialdata.allchoices = repmat({1:1:180},nTrials,1);

end
