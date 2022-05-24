% Script for making an axis label demonstrating how the difficulty might
% increase as the distance of the closest distractor becomes further away

clear, clc, clear

cleandata = analysis('dirname', {'210905_112706_Pollux'},...
    'stimhistfig',1,'makefigs',0); % arbitrary file, just to quickly get colvals

stimCols = cleandata.trialdata.stimCols{1,1};
pltCols = double(LuvTosRGB(stimCols))/255;

%%
figure, hold on

cueCol = 28;
pltSize = 100; % FAO Audrey

scatter(-1,2,...
    pltSize, pltCols(cueCol,:), "filled")
scatter(1:31, ones(31,1)*2,...
    pltSize, pltCols(cueCol+1:cueCol+31,:), "filled")

scatter(-1,1,...
    pltSize, pltCols(cueCol,:), "filled")
scatter(1:31, ones(31,1), pltSize,...
    pltCols([cueCol-1:-1:1,64:-1:61],:), "filled")

axis equal
