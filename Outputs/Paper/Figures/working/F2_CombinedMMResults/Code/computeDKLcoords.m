clc, clear, close all

% Script for computing the DKL axis locations for F2C
% Modified/extended copy of CausalGlobs\protocol\calibration\20230108\Analysis20230108.m

%%

load(['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,...
    'Data',filesep,'stimuliMeasurements',filesep,'StimuliMeasurements20230108.mat'])

%%

%load gammaCorrectionValues_18-Apr-2021.mat
%whiteXYZ = xyYToXYZ(LumValues.white(end).xyYcie');
whiteXYZ = XYZ(:,430)*2; % where the measurement had looped and I was measuring gray. And then times 2 to roughly get us to the white point. Need to read in the actual white point that is specified. !!!!!!!!!!!!!!!!!!!!!!!

% load C:\Users\cege-user\Documents\Kofiko\calibrationData\conversionMatrices.mat

Luv = XYZToLuv(XYZ,whiteXYZ);

xyY = XYZToxyY(XYZ);
%%
scatter3(Luv(2,:),Luv(3,:),1:size(Luv,2))

%scatter3(Luv(2,:),Luv(3,:),Luv(1,:))

axis auto%equal

%% Independently compute XYZ from spectra

S_SPD = [380,4,101];

load T_xyz1931.mat T_xyz1931 S_xyz1931
T_xyz1931_interp = SplineCmf(S_xyz1931,T_xyz1931,S_SPD);

figure, hold on
plot(SToWls(S_xyz1931),T_xyz1931,'--','LineWidth',2)
plot(SToWls(S_SPD),T_xyz1931_interp,'--','LineWidth',3)
legend

XYZ_indep = T_xyz1931_interp*SPD;
% XYZ = XYZ/max(XYZ(2,:));


for i = 1:3
    figure,
    scatter3(XYZ(i,:),XYZ_indep(i,:),1:size(XYZ,2),...
        'k','filled','MarkerFaceAlpha',0.2)
    view(2)
end
%% Chroma

C = sqrt(Luv(2,:).^2 + Luv(3,:).^2);

figure, 
plot(C)

%% Cut out gray section in the middle

figure,
% plot(SPD,'k')
plot(SPD(:,550:600),'k')
SPDbackground = mean(SPD(:,550:600),2);

SPD = SPD(:,[1:350,700:end]);
figure,
% plot(SPD,'k')
plot(SPD(:,1:51),'k')

XYZ_background = mean(XYZ_indep(:,550:600),2);
XYZ_indep = XYZ_indep(:,[1:350,700:end]);

xyY_background = XYZToxyY(XYZ_background);

figure,
scatter3(Luv(2,:),Luv(3,:),Luv(1,:))
axis square
view(2)

Luv = Luv(:,[1:350,700:end]);
figure,
scatter3(Luv(2,:),Luv(3,:),Luv(1,:))
axis square
view(2)

%%

% stimCols = generateStimCols('nBig',64);
addpath(genpath(['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep]))
stimCols = generateStimCols('nBig',64,'showFig',true);
hold on
scatter(stimCols(1,:),stimCols(2,:))

stimColsLuv = [ones(1,size(stimCols,2))*76;stimCols];

%% Compare to specified colors

rng(42) % fix random number generator for reproducibility
[k_IDX, C, SUMD, D] = kmeans(Luv(2:3,:)', 64, 'Display','iter','Replicates',10);


figure, hold on
for i = 1:64
    plot(Luv(2,k_IDX==i),Luv(3,k_IDX==i),'LineWidth',3)

    s(i,:) = std(Luv(:,k_IDX==i),[],2);
    text(Luv(2,find(k_IDX==i, 1)),Luv(3,find(k_IDX==i, 1)),num2str(i))
    
end

hold on
scatter(stimCols(1,:),stimCols(2,:))
% for i = 1:64
%     text(stimCols(1,i),stimCols(2,i),num2str(i))
% end

axis equal

% for i = 1:64
%     [~,minInd(i)] = min(mean(abs(stimCols(:,i) - C')));
%     comet([stimCols(1,i),C(minInd(i),1)],[stimCols(2,i),C(minInd(i),2)])
%     C(minInd(i),:) = Inf; % hacky method that means it is essentially sampling without replacement
% end

% figure, hold on
% for i = 2
%     for j = 1:size(SPD,2)
%         if k_IDX(j) == i
%             plot(SPD(:,j))
%             disp(k_IDX(j))
%         end
%     end
% end

%% order by angle

% Compute angles for each k-means group
angs = rad2deg(cart2pol(C(:,1),C(:,2)));
angs(angs<-2) = angs(angs<-2) + 360; % This would normally be 0 rather than -2, this is just hacky since I can see that the first value is roughly -1.8

% sort by angle, and store the order of the remapping
[angs_sorted,angs_sorted_index] = sort(angs);

stimulusIndex = NaN(1,size(SPD,2));
for i = 1:length(stimulusIndex)
    stimulusIndex(i) = find(k_IDX(i) == angs_sorted_index);
end

% % calculate the index for each recording, in terms of which stimulus it corresponds to
% ind2 = angs_sorted_index(k_IDX); 


%%
% for j = 1:64 % for each stimulus
%     figure, hold on
%     for i = 1:size(SPD,2) % go through every SPD        
%         if stimulusIndex(i) == j % if the index for that recording matches the requested stimulus
%             plot(SPD(:,i)) % plot it
%         end
%     end
% end

figure, plot(SPD)

for j = 1:64
    SPD_avs(:,j) = mean(SPD(:,stimulusIndex == j),2);
end

figure, plot(SPD_avs)

%% Compute DKL values

load T_cones_sp.mat T_cones_sp S_cones_sp
load T_xyzJuddVos.mat T_xyzJuddVos S_xyzJuddVos
T_Y = 683*T_xyzJuddVos(2,:);
T_Y = SplineCmf(S_xyzJuddVos,T_Y,S_cones_sp);

bgLMS   = T_cones_sp * SplineSpd(S_SPD,SPDbackground,S_cones_sp);
% LMS     = (T_cones_sp * SplineSpd(S_SPD,SPD,S_cones_sp)) - bgLMS;
LMSinc  = (T_cones_sp * SplineSpd(S_SPD,SPD_avs,S_cones_sp)) - bgLMS;

[M_ConeIncToDKL,LMLumWeights] = ComputeDKL_M(bgLMS,T_cones_sp,T_Y);

DKL = M_ConeIncToDKL*LMSinc;

% sRGBlin = XYZToSRGBPrimary(XYZ_indep/(XYZ_background(2)*2)); % TODO Hacky whitepoint, visualization only though
% sRGB = uint8(SRGBGammaCorrect(sRGBlin,0)');

% figure,
% scatter3(Luv(2,:),Luv(3,:),Luv(1,:),...
%     [],double(sRGB)/255,'filled','MarkerEdgeColor','k')
% axis square
% view(2)

figure,
scatter3(DKL(2,:),DKL(3,:),DKL(1,:))
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')
title('DKL')

figure,
scatter3(DKL(2,:),DKL(3,:),1:64)
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('StimIndex')
title('DKL')

%% Compute locations of axes

lm_temp = DKL(3,:);
lm_temp(:,DKL(2,:) < 0) = inf;
[~,lm_plus]  = min(abs(lm_temp));
lm_temp = DKL(3,:); 
lm_temp(:,DKL(2,:) > 0) = inf;
[~,lm_minus]  = min(abs(lm_temp));

s_temp = DKL(2,:);
s_temp(:,DKL(3,:) < 0) = inf;
[~,s_plus]  = min(abs(s_temp));
s_temp = DKL(2,:); 
s_temp(:,DKL(3,:) > 0) = inf;
[~,s_minus]  = min(abs(s_temp));

% reassurance
figure, hold on
scatter3(DKL(2,lm_plus),DKL(3,lm_plus),lm_plus,'r','filled','DisplayName','lm\_plus')
scatter3(DKL(2,lm_minus),DKL(3,lm_minus),lm_minus,'g','filled','DisplayName','lm\_minus')
scatter3(DKL(2,s_plus),DKL(3,s_plus),s_plus,'b','filled','DisplayName','s\_plus')
scatter3(DKL(2,s_minus),DKL(3,s_minus),s_minus,'y','filled','DisplayName','s\_minus')
legend('AutoUpdate','off')
scatter3(DKL(2,:),DKL(3,:),1:64,'k')

xlabel('L-M')
ylabel('S-(L+M)')
zlabel('StimIndex')
title('DKL')

xline(0)
yline(0)


saveas(gcf,['../','stimInDKL_',datestr(now,'yymmdd-HHMMSS'),'.svg'])


disp(lm_plus)
disp(lm_minus)
disp(s_plus)
disp(s_minus)


stimCols = generateStimCols('nBig',64,'showFig',true);
s1 = scatter(stimCols(1,lm_plus),stimCols(2,lm_plus),'r','filled','DisplayName','lm\_plus');
s2 = scatter(stimCols(1,lm_minus),stimCols(2,lm_minus),'g','filled','DisplayName','lm\_minus');
s3 = scatter(stimCols(1,s_plus),stimCols(2,s_plus),'b','filled','DisplayName','s\_plus');
s4 = scatter(stimCols(1,s_minus),stimCols(2,s_minus),'y','filled','DisplayName','s\_minus');
legend([s1,s2,s3,s4],'location','best')
title('CIELUV')

saveas(gcf,['../','stimInLUVwithDKLpolesHighlighted_',datestr(now,'yymmdd-HHMMSS'),'.svg'])

%% Parallel conversion: using XYZtoLMS conversion matrix instead of independent spectral measurements

XYZtoLMS_HPEee = [0.38971, 0.68898, -0.07868;...
                 -0.22981, 1.18340,  0.04641;...
                  0,       0,        1]; 
%Hunt-Pointer-Estevez (equi-energy) - https://en.wikipedia.org/wiki/LMS_color_space#Hunt,_RLAB

bgLMS2  = XYZtoLMS_HPEee*XYZ_background;
LMSinc2  = XYZtoLMS_HPEee*XYZ;

[M_ConeIncToDKL2] = ComputeDKL_M(bgLMS2,T_cones_sp,T_Y);

DKL2 = M_ConeIncToDKL2*LMSinc2;

figure,
scatter3(DKL2(2,:),DKL2(3,:),DKL2(1,:))
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')
title('DKL (by HPE)')

%%

load('T_cones_ss2')
load('T_ss2000_Y2.mat')
T_ss2000_Y2 = 683*T_ss2000_Y2;


XYZtoLMS_SS = [ 0.210576,   0.855098,  -0.0396983;
               -0.417076,   1.177260,   0.0786283;
                0,          0,          0.5168350];
% Stockman Sharpe 2000 (https://en.wikipedia.org/wiki/LMS_color_space#Stockman_&_Sharpe_(2000))


bgLMS3  = XYZtoLMS_SS*XYZ_background;
LMSinc3 = XYZtoLMS_SS*XYZ;

[M_ConeIncToDKL3] = ComputeDKL_M(bgLMS3,T_cones_ss2,T_ss2000_Y2);

DKL3 = M_ConeIncToDKL3*LMSinc3;

figure,
scatter3(DKL3(2,:),DKL3(3,:),DKL3(1,:))
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')
title('DKL (by SS)')


%% Compute background luminance

lumFunc = 683*T_xyz1931(2,:);

SPDbackground_int = SplineSpd(S_SPD,SPDbackground,S_xyz1931);
figure, hold on
plot(SToWls(S_SPD),SPDbackground)
plot(SToWls(S_xyz1931),SPDbackground_int)

lumBackground = lumFunc*SPDbackground_int;









