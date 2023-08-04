% Draw chromaticity diagram

% Requires PTB
% Requires TCC_AFC

clear, clc, close all

addpath(genpath(['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Analyses']))

%%

figure

% DrawChromaticity('upvp')

% Can I add u*v* as second axes? !!!!
% Should I add gamut? !!!!

%% Add stimuli

Lstar = 76.0693;
whitepoint = [0.3,0.4]; % Need to grab the actual white point !!!!

stimCols = [generateStimCols('nBig',64),[0;0]];
% stimCols_upvp = stimCols./(13*Lstar) + whitepoint';
stimCols_sRGB = im2double(LuvTosRGB([repelem(Lstar, 65); stimCols(1,:); stimCols(2,:)]))';

hold on
scatter(stimCols(1,:),stimCols(2,:),100,stimCols_sRGB',"filled")

xticks([-40,0,40])
yticks([-40,0,40])

xlabel('u*')
ylabel('v*')

%%

axis equal square

saveas(gcf,fullfile(['../','chromaticity', datestr(now,'yymmdd'), '.svg']))

