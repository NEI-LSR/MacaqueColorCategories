clear, clc, close all

%% Load similarity matrix

repoHomeDir = ['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..'];
addpath(genpath([repoHomeDir,filesep,'Analyses']))
% ModelDir = [repoHomeDir,filesep,'Analyses',filesep,'TCCModels',filesep,'Combined',filesep];

filename = ['..',filesep,'..',filesep,'..',filesep,'F4_TCCResults',filesep,'Code',filesep,'combined_TCC-FreeSimilarityMatrix_231116-190823']; % TODO Move this file
load([filename,'.mat'],'x')

plotSimilarityMatrix(x);

nBig = sqrt(length(x));
sm = reshape(x,[nBig, nBig]);

%% Extract columns

% which way round is this matrix?
% sm(i,:) is a column of the similarity matrix, aka a cue

% circular rotate
for i = 1:nBig
    sm_cs(i,:) = circshift(sm(i,:),(nBig - i) + (nBig/2));
    % sm_cs(i,:) = circshift(sm(i,:), nBig - i);
end
figure, 
imagesc(sm_cs');
colorbar
axis image
colormap('gray')

%% Compute x-axis values

interval = 360/nBig;
hueAngles = -180+interval:interval:180;

[cart,pol] = generateStimCols('nBig',nBig);
for i = 1:nBig
    dE(i) = sqrt(... % deltaE, aka euclidian distance
        ((cart(1,1)-cart(1,i)).^2)...
        +...
        ((cart(2,1)-cart(2,i)).^2));
end
dE_cs = circshift(dE,31);

%%
f = figure; hold on
yyaxis left
plot(1:nBig,pol(1,:))

xlabel('Stimulus Index')
ylabel('Angular distance (degrees)')

yyaxis right
plot(1:nBig,dE);

xlim([1 33])
ylabel('Euclidian Distance (Delta E)')

%%

saveas(f,['../','distanceMetricComparisonA','.svg'])

%% Plot individual columns

stimCols = generateStimCols('nBig',nBig);
stimCols_sRGB = LuvTosRGB([ones(1,nBig)*76.0693;stimCols]);

figure, hold on
for i = 1:nBig
    plot(hueAngles,sm_cs(i,:),...
        'Color',[stimCols_sRGB(i,:),255/3],...
        'HandleVisibility','off')
end

xlim([min(hueAngles),max(hueAngles)])

xlabel('Hue angle (degrees)')
ylabel('Similarity')

%% Compute and plot average

plot(hueAngles,mean(sm_cs,1),'kx-','DisplayName','Mean')
plot(hueAngles,median(sm_cs,1),'ko-','DisplayName','Median')
legend

xline(0,'HandleVisibility','off')

%% Fit gaussians

gaussEqn = @(sd,x) exp(-((x.^2)/(2*sd^2)));
% gaussEqn = 'exp(-((x^2)/(2*sd^2)))';

% Hue angle
x0 = 30;
f_hueAngle = fit(hueAngles', median(sm_cs,1)', ...
    gaussEqn,'start',x0);

y_hueAngle = feval(f_hueAngle,hueAngles);
plot(hueAngles,y_hueAngle,'k:','LineWidth',2,...
    'DisplayName','Gaussian: Hue Angle')

% Cartesian

x0 = 30;
f_dE = fit(dE_cs', median(sm_cs,1)', ...
    gaussEqn,'start',x0);

y_dE = feval(f_dE,dE_cs);
plot(hueAngles,y_dE,'k--','LineWidth',2,...
    'DisplayName','Gaussian: Delta E')

%%
saveas(gcf,['../','distanceMetricComparisonB','.svg'])