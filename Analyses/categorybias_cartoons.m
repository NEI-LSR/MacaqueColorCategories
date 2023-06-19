%clear, clc, close all


%%
bigN = 64;
smallN = 4;
hue_angle = 0:360/bigN:360;
hue_angle = hue_angle(1:bigN);

all_cues = (1:bigN)';
degrees = 360/length(all_cues);

% Colors
stimCols = generateStimCols('nBig',bigN,'sat',37);
rotVal = 360/(bigN*2);
rotationMatrix = [cosd(rotVal),-sind(rotVal);sind(rotVal),cosd(rotVal)];
stimCols_biased = rotationMatrix * stimCols;
stimCols_sRGB = LuvTosRGB([repelem(76.0693, bigN); stimCols]);
rstimCols_sRGB = LuvTosRGB([repelem(76.0693, bigN); stimCols_biased]);
colvals = im2double(stimCols_sRGB);

%% Pie chart figure
hue_angle = 0:360/bigN:360;
hue_angle = hue_angle(1:bigN);

figure('WindowState', 'maximized');
hold on
pie(repelem(degrees,bigN));% pie chart w/ bigN equally sized slices
colormap(ones(size(rstimCols_sRGB)));
ax = gca;
delete(ax.Children([1, 1:2:bigN*2])) % stop displaying % for each slice
for i = 1:bigN
    ax.Children(i).EdgeAlpha = 0; % get rid of lines between slices
end
for i = 1:bigN
    ax.Children(i).FaceAlpha = 1;
end
ax.View = [90 90];


rotated_colvals = im2double(rstimCols_sRGB);

shift_colvals = [colvals(48:end,:); colvals(1:47,:)];
[cart,~] = generateStimCols('nBig',64); % generate values to plot cues

scatter(cart(1,:)./36,cart(2,:)./36,185,shift_colvals,'filled'); % plot all cues around pie chart

set(gca,'visible','off')
axis equal

rad_angle = deg2rad(hue_angle); %hue angles in radians
axes3 = axes('Position', [0.318 0.1 0.4 0.835]);

%%

% % none
% polarplot(rad_angle([1:end 1]), zeros(length(rad_angle([1:end 1])))+25,'k','LineWidth',5);

% % human
% 
% ngps = [32,146,235,326]; %negcrosspoints
% prec = 1000; %precisions 
% scaler = 10; % to match bias
% 
% x_all = [linspace(-(360-ngps(end)),ngps(1),prec),...
%     linspace(ngps(1),ngps(2),prec),...
%     linspace(ngps(2),ngps(3),prec),...
%     linspace(ngps(3),ngps(4),prec),...
%     linspace(ngps(4),360+ngps(1),prec)];
% 
% y = -sin(linspace(0,1,prec)*2*pi);
% 
% y_all = [y,y,y,y,y] * scaler;
% 
% plot(x_all,y_all)
% xlim([0 360])
% 
% polarplot(deg2rad(x_all),y_all +25,'k','LineWidth',5);
% hold on
% for k = 1:length(ngps)
%     polarplot([deg2rad(ngps(k)) deg2rad(ngps(k))],[0 100],'Color','k');
% end

% % unique hues
% 
% ngps = [30,100,150,250]; %negcrosspoints
% prec = 1000; %precisions 
% scaler = 10; % to match bias
% 
% x_all = [linspace(-(360-ngps(end)),ngps(1),prec),...
%     linspace(ngps(1),ngps(2),prec),...
%     linspace(ngps(2),ngps(3),prec),...
%     linspace(ngps(3),ngps(4),prec),...
%     linspace(ngps(4),360+ngps(1),prec)];
% 
% y = -sin(linspace(0,1,prec)*2*pi);
% 
% y_all = [y,y,y,y,y] * scaler;
% 
% plot(x_all,y_all)
% xlim([0 360])
% 
% polarplot(deg2rad(x_all),y_all +25,'k','LineWidth',5);
% hold on
% for k = 1:length(ngps)
%     polarplot([deg2rad(ngps(k)) deg2rad(ngps(k))],[0 100],'Color','k');
% end

% % dkl
% 
% ngps = [0,78.75,191.25,258.75]; %negcrosspoints
% prec = 1000; %precisions 
% scaler = 10; % to match bias
% 
% x_all = [linspace(-(360-ngps(end)),ngps(1),prec),...
%     linspace(ngps(1),ngps(2),prec),...
%     linspace(ngps(2),ngps(3),prec),...
%     linspace(ngps(3),ngps(4),prec),...
%     linspace(ngps(4),360+ngps(1),prec)];
% 
% y = -sin(linspace(0,1,prec)*2*pi);
% 
% y_all = [y,y,y,y,y] * scaler;
% 
% plot(x_all,y_all)
% xlim([0 360])
% 
% polarplot(deg2rad(x_all),y_all +25,'k','LineWidth',5);
% 
% hold on
% for k = 1:length(ngps)
%     polarplot([deg2rad(ngps(k)) deg2rad(ngps(k))],[0 100],'Color','k');
% end

%%
rlim([0 50]);
hold on
polarplot(rad_angle([1:end 1]), zeros(length(rad_angle)+1,1)+25,'LineStyle','--','Color','k');
ax = gca;
ax.Color = 'none';
ax.ThetaTickLabel = {'0','','','90','','','180','','','270','',''};
ax.RTick = 0:5:50;
rticklabels({'','-20{\circ}','','-10{\circ}','','0{\circ}','','10{\circ}','','20{\circ}',''});
ax.RAxisLocation = 235;
ax.FontSize = 8;
ax.Units = 'normalized';
ax.Position = [0.3734375,0.212952799121844,0.2890625,0.609220636663008];



