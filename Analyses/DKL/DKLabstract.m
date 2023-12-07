% LabToDKL

clear, clc, close all

% Lab = [60,60,60,60,60;...
%     0,20,0,-20,0;...
%     0,0,20,0,-20];

addpath(genpath(['..',filesep]))
[cart,pol] = generateStimCols('nBig',64,'sat',20);
LXX = [ones(1,64)*60;cart];

LXX(:,end+1) = [60;0;0];

% figure,
% scatter(Lab(2,:),Lab(3,:),Lab(1,:))

whiteXYZ = [95.04;100;108.88];
% XYZ = LabToXYZ(LXX,whiteXYZ);
XYZ = LuvToXYZ(LXX,whiteXYZ);

load T_cones_sp
load T_xyz1931.mat
S_cones = S_cones_sp;
T_cones = T_cones_sp;
T_xyz = T_xyz1931;
S_xyz = S_xyz1931;

T_Y = 683*T_xyz(2,:);
S_Y = S_xyz;
T_Y = SplineCmf(S_Y,T_Y,S_cones);

T_xyz = SplineCmf(S_xyz,T_xyz,S_cones);

M_XYZToLMS = (T_xyz'\T_cones')';  
T_cones_chk = M_XYZToLMS*T_xyz;

figure, hold on
plot(T_cones','k')
plot(T_cones_chk','r--')

bgLMS  = M_XYZToLMS * XYZ(:,end);
LMSinc = M_XYZToLMS * (XYZ - XYZ(:,end));

[M_ConeIncToDKL] = ComputeDKL_M(bgLMS, T_cones, T_Y);
DKL = M_ConeIncToDKL * LMSinc;

figure,
scatter3(DKL(2,:),DKL(3,:),DKL(1,:),'DisplayName','DKL (self-computed)')
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')
view(2)

legend

%%

lm_temp = DKL(3,1:64);
lm_temp(:,DKL(2,:) < 0) = inf;
[~,lm_plus]  = min(abs(lm_temp));
lm_temp = DKL(3,1:64); 
lm_temp(:,DKL(2,:) > 0) = inf;
[~,lm_minus]  = min(abs(lm_temp));

s_temp = DKL(2,1:64);
s_temp(:,DKL(3,:) < 0) = inf;
[~,s_plus]  = min(abs(s_temp));
s_temp = DKL(2,1:64); 
s_temp(:,DKL(3,:) > 0) = inf;
[~,s_minus]  = min(abs(s_temp));

figure, hold on
scatter3(DKL(2,lm_plus),DKL(3,lm_plus),lm_plus,'r','filled','DisplayName','lm\_plus')
scatter3(DKL(2,lm_minus),DKL(3,lm_minus),lm_minus,'g','filled','DisplayName','lm\_minus')
scatter3(DKL(2,s_plus),DKL(3,s_plus),s_plus,'b','filled','DisplayName','s\_plus')
scatter3(DKL(2,s_minus),DKL(3,s_minus),s_minus,'y','filled','DisplayName','s\_minus')
legend('AutoUpdate','off')
scatter3(DKL(2,:),DKL(3,:),1:65,'k')

xlabel('L-M')
ylabel('S-(L+M)')
zlabel('StimIndex')
title('DKL')

xline(0)
yline(0)

disp(lm_plus)
disp(lm_minus)
disp(s_plus)
disp(s_minus)
