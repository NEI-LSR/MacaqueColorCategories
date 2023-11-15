clear, clc, close all

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))

%%

filename{1} = '210422--211012_Pollux_data';
filename{2} = '210517--211108_Castor_data';
filename{3} = '210428--210609_Buster_data';
filename{4} = '220322--220823_Morty_data';

% TODO Replace with CSVs

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses'))

for participant = 1:length(filename)

    load(['C:\Users\cege-user\Documents\MacaqueColorCategories\Data\',...
        filename{participant},'.mat'])
    cleandata.trialdata = trialdata;

    difficulty_psychometric(cleandata,filename{participant},true); % all on one graph
    % this saves each intermediary figure too
    % so they need manually deleting afterwards

    difficulty_psychometric(cleandata,filename{participant});

end

