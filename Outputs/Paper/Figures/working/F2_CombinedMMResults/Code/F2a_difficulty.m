clear, clc, close all

% before running, `cd` to where this file is

filename{1} = '210422--211012_Pollux_data';
filename{2} = '210517--211108_Castor_data';
filename{3} = '210428--210609_Buster_data';
filename{4} = '220322--220823_Morty_data';

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
addpath(genpath([repoHomeDir,'Analyses']))


%%

for participant = 1:length(filename)

    load([repoHomeDir,filesep,'Data',filesep,filename{participant},'.mat'],'cleandata');
    cleandata.nBig = 64;

    difficulty_psychometric(cleandata,filename{participant},true); % all on one graph
    % this saves each intermediary figure too
    % so they need manually deleting afterwards

    difficulty_psychometric(cleandata,filename{participant});

end

%% to get the number of cues per animal that have been subsampled to match Buster

% cues = cell2mat(cleandata.trialdata.cues);
% 
% for i = 1:64
%     t(i,1) = sum(cues(1:24526) == i);
% end
% 
% for i = 1:64
%     t(i,2) = sum(cues(24527:49052) == i);
% end
% 
% for i = 1:64
%     t(i,3) = sum(cues(49053:73578) == i);
% end
% 
% for i = 1:64
%     t(i,4) = sum(cues(73579:end) == i);
% end
% 
% disp(sum(t(:)) == size(cues,1)) % sanity check
% 
% 

