clear, clc, close all

filename{1} = '210422--211012_Pollux_data';
filename{2} = '210517--211108_Castor_data';
filename{3} = '210428--210609_Buster_data';
filename{4} = '220322--220823_Morty_data';

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath([repoHomeDir,'Analyses']))


%%

% TODO Add in CSV option

for participant = 1:length(filename)

    cleandata = load([repoHomeDir,filesep,'Data',filesep,filename{participant},'.mat']);

    difficulty_psychometric(cleandata,filename{participant},true); % all on one graph
    % this saves each intermediary figure too
    % so they need manually deleting afterwards

    difficulty_psychometric(cleandata,filename{participant});

end

