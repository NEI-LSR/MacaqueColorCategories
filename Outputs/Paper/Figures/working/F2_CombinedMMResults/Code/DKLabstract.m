% LabToDKL

clear, clc, close all

% Lab = [60,60,60,60,60;...
%     0,20,0,-20,0;...
%     0,0,20,0,-20];

addpath(genpath(['..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'..',filesep,'Analyses']))
Lab = [ones(1,64)*60;generateStimCols('nBig',64,'sat',20)];

Lab(:,end+1) = [60;0;0];

% figure,
% scatter(Lab(2,:),Lab(3,:),Lab(1,:))

whiteXYZ = [95.04;100;108.88];
XYZ = LabToXYZ(Lab,whiteXYZ);

load T_cones_sp
load T_xyzJuddVos
S_cones = S_cones_sp;
T_cones = T_cones_sp;
T_Y = 683*T_xyzJuddVos(2,:);
S_Y = S_xyzJuddVos;
T_Y = SplineCmf(S_Y,T_Y,S_cones);

XYZtoLMS_HPEee = [0.38971, 0.68898, -0.07868;...
                 -0.22981, 1.18340,  0.04641;...
                  0,       0,        1]; 
%Hunt-Pointer-Estevez (equi-energy) - https://en.wikipedia.org/wiki/LMS_color_space#Hunt,_RLAB

M_XYZToLMS = (T_xyzJuddVos'\T_cones')';  
T_cones_chk = M_XYZToLMS*T_xyzJuddVos;

% figure, hold on
% plot(T_cones','k')
% plot(T_cones_chk','r--')

bgLMS_HPEee  = XYZtoLMS_HPEee * XYZ(:,1);
LMSinc_HPEee = XYZtoLMS_HPEee * (XYZ - XYZ(:,1));
bgLMS_sc  = M_XYZToLMS * XYZ(:,end);
LMSinc_sc = M_XYZToLMS * (XYZ - XYZ(:,end));

[M_ConeIncToDKL_HPEee] = ComputeDKL_M(bgLMS_HPEee, T_cones, T_Y);
DKL_HPEee = M_ConeIncToDKL_HPEee * LMSinc_HPEee;
[M_ConeIncToDKL_sc] = ComputeDKL_M(bgLMS_sc, T_cones, T_Y);
DKL_sc = M_ConeIncToDKL_sc * LMSinc_sc;

figure, hold on
scatter3(DKL_HPEee(2,:),DKL_HPEee(3,:),DKL_HPEee(1,:),'DisplayName','DKL (via HPEee)')
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')

scatter3(DKL_sc(2,:),DKL_sc(3,:),DKL_sc(1,:),'DisplayName','DKL (self-computed)')
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M')
view(2)

legend
