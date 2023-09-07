clear, clc, close all

p = {'Castor','Pollux','Buster','Morty'};

for i = 1:length(p)
    plotbar_NLL_AIC_BIC(p{i},'../')
    close all
end