% Draw chromaticity diagram

% Requires PTB
% Requires TCC_AFC

clear, clc, close all

%%

DrawChromaticity('upvp')
% Can I add u*v* as second axes? !!!!
% Should I add gamut? !!!!

%% Add stimuli

Lstar = 76.0693;
whitepoint = [0.3,0.4]; % Need to grab the actual white point !!!!

stimCols = generateStimCols('nBig',64);
stimCols_upvp = stimCols./(13*Lstar) + whitepoint';
stimCols_sRGB = im2double(LuvTosRGB([repelem(Lstar, 64); stimCols(1,:); stimCols(2,:)]))';

hold on
scatter(stimCols_upvp(1,:),stimCols_upvp(2,:),[],stimCols_sRGB',"filled")

xlim([0, 0.6])
ylim([0, 0.6])

%%

f = gcf;

f.Position = [100 100 175 175];

formatFig(gcf, [3 3], 'nature') % bc6:\MakeFigure
saveas(gcf,fullfile(['chromaticity', datestr(now,'yymmdd'), '.svg']))

