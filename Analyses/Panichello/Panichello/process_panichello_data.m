% Load Panichello data
filepath = 'C:\Users\selwynhm\Documents\MATLAB\nAFC\analysis\Panichello\';
filename = 'human.mat';

load([filepath,filename])

cues = cleandata.trialdata.cues;
chosen = cleandata.trialdata.chosen;

chosen(chosen > 360) = chosen(chosen > 360) - 360; % for human data only

choices = repmat({1:1:360},length(cues),1);

cleandata.trialdata.cues = num2cell(cues);
cleandata.trialdata.chosen = num2cell(chosen);
cleandata.trialdata.choices = choices;
cleandata.trialdata.dirname = 'Panichello_Humans';
cleandata.trialdata.paradigm = 'Panichello_Humans';
cleandata.trialdata.allchoices = num2cell(1:360);
