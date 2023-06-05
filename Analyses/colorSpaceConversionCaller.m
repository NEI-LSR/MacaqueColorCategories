% For generating n points that are uniformally sampled in the behaviorally
% derived colorspace.

clc, clear, close all

n = 1000; % number of points

%%
requestedHueAngles = linspace(0,360,n+1);
requestedHueAngles = requestedHueAngles(1:end-1);
indexedLocation = colorSpaceConversion(requestedHueAngles);

%%

[stimColsCIELUVcart,stimColsCIELUVpolar] = generateStimCols('nBig',64);
interval = stimColsCIELUVpolar(1,2)-stimColsCIELUVpolar(1,1);
for sc = 1:n
    stimColsMACBEHpolar(sc) = stimColsCIELUVpolar(1,indexedLocation(sc,1)) + (interval*indexedLocation(sc,2));
end

%% Convert to cartesian chromaticities for plotting

[stimCols_MACBEHcart(1,:),stimCols_MACBEHcart(2,:)] = pol2cart(deg2rad(stimColsMACBEHpolar), stimColsCIELUVpolar(2,1));

stimCols_sRGB = LuvTosRGB([repelem(76.0693, n); stimCols_MACBEHcart(1,:); stimCols_MACBEHcart(2,:)]);
colvals = im2double(stimCols_sRGB);


%%

[stimColsPlotting,stimColsPlottingPolar] = generateStimCols('nBig',n,'showFig',true);
axis off
title('CIELUV')

figure,
hold on
axis equal
axis off
title('Behav. Defined')

scatter(stimColsPlotting(1,:),stimColsPlotting(2,:),[],colvals,"filled")

%% 

figure, 
hold on
plot(stimColsPlottingPolar(1,:),'DisplayName','CIELUV')
plot(stimColsMACBEHpolar,'DisplayName','Behav. Defined')
xlabel('n_{sc}')
ylabel('hue angle')
legend('Location','best')
