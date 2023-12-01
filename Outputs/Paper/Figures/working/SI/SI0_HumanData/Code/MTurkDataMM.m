% Conway Lab MTurk Data

clc, clear, close all

%%

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep];
dataDir = [repoHomeDir,'Data'];

addpath(genpath(repoHomeDir));

%%

load([dataDir,filesep,'211112--220507_MTurk.mat'])

%%

Lab = 0;
lengthOfSlidingWindow = [];
includeCorrect = false;

model = fitMixtureModel(cleandata,Lab,lengthOfSlidingWindow,includeCorrect);

%%

whichFigures.MixMod_polar    = true;
whichFigures.MixMod_linear   = true;

plotMixtureModel(model,...
    whichFigures,'MTurk')


