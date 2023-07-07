%% category comparison script

% category_data = 'AllData/category_data.mat';
% 
 %load(category_data);
 
 load category_data.mat

offset = 0;

figure;
hold on

% Pollux
polarAngs = (category_data.Pollux.centers)'; %Polar Angles
[a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
stimCols = [a;b];
stimCols_sRGB = LuvTosRGB([repelem(76.0693, length(polarAngs)); stimCols]);
colvals = im2double(stimCols_sRGB);

p3_err_neg = category_data.Pollux.centers - category_data.Pollux.confidence(1:2:end);
p3_err_pos = mod(category_data.Pollux.confidence(2:2:end) - category_data.Pollux.centers, 360);
for i = 1:length(p3_err_neg)
errorbar(category_data.Pollux.centers(i)- offset, repelem(7, 1),p3_err_neg(i), p3_err_pos(i), 'horizontal','ko','MarkerFaceColor',colvals(i,:),'MarkerEdgeColor',colvals(i,:));
end

% Castor
polarAngs = (category_data.Castor.centers)'; %Polar Angles
[a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
stimCols = [a;b];
stimCols_sRGB = LuvTosRGB([repelem(76.0693, length(polarAngs)); stimCols]);
colvals = im2double(stimCols_sRGB);

p4_err_neg = category_data.Castor.centers - category_data.Castor.confidence(1:2:end);
p4_err_pos = mod(category_data.Castor.confidence(2:2:end) - category_data.Castor.centers, 360);
for i = 1:length(p4_err_neg)
errorbar(category_data.Castor.centers(i)- offset, repelem(6,1),p4_err_neg(i), p4_err_pos(i), 'horizontal','ko', 'MarkerEdgeColor',colvals(i,:),'MarkerFaceColor',colvals(i,:));
end
% Buster
polarAngs = (category_data.Buster.centers)'; %Polar Angles
[a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
stimCols = [a;b];
stimCols_sRGB = LuvTosRGB([repelem(76.0693, length(polarAngs)); stimCols]);
colvals = im2double(stimCols_sRGB);

p5_err_neg = category_data.Buster.centers - category_data.Buster.confidence(1:2:end);
p5_err_pos = mod(category_data.Buster.confidence(2:2:end) - category_data.Buster.centers, 360);
for i = 1:length(p5_err_neg)
errorbar(category_data.Buster.centers(i)- offset, repelem(5,1),p5_err_neg(i), p5_err_pos(i), 'horizontal', 'ko', 'MarkerEdgeColor',colvals(i,:),'MarkerFaceColor',colvals(i,:));
end
errorbar(-category_data.Buster.centers(2) + offset, repelem(5,1),p5_err_neg(2), p5_err_pos(2), 'horizontal', 'k', 'MarkerEdgeColor',colvals(2,:),'MarkerFaceColor',colvals(2,:));

% Morty - 3 centers, 2 interval edges? 
avg_center = (category_data.Morty.centers(2) + category_data.Morty.centers(3))/2;
category_data.Morty.centers = [category_data.Morty.centers(1); avg_center];
polarAngs = (category_data.Morty.centers)'; %Polar Angles
[a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
stimCols = [a;b];
stimCols_sRGB = LuvTosRGB([repelem(76.0693, length(polarAngs)); stimCols]);
colvals = im2double(stimCols_sRGB);
p6_err_neg = category_data.Morty.centers - category_data.Morty.confidence(1:2:end);
p6_err_pos = category_data.Morty.confidence(2:2:end) - category_data.Morty.centers;
for i = 1:length(p6_err_neg)
errorbar(category_data.Morty.centers(i) - offset, repelem(4,1),p6_err_neg(i), p6_err_pos(i), 'horizontal', 'ko', 'MarkerEdgeColor',colvals(i,:),'MarkerFaceColor',colvals(i,:));
end

% MechTurk
polarAngs = (category_data.MechTurk.centers)'; %Polar Angles
[a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
stimCols = [a;b];
stimCols_sRGB = LuvTosRGB([repelem(76.0693, length(polarAngs)); stimCols]);
colvals = im2double(stimCols_sRGB);

p2_err_neg = category_data.MechTurk.centers - category_data.MechTurk.confidence(1:2:end);
p2_err_pos = mod(category_data.MechTurk.confidence(2:2:end) - category_data.MechTurk.centers, 360);
for i = 1:length(p2_err_neg)
errorbar(category_data.MechTurk.centers(i)- offset, repelem(3,1),p2_err_neg(i), p2_err_pos(i), 'horizontal', 'ko','MarkerEdgeColor',colvals(i,:),'MarkerFaceColor',colvals(i,:));
end

% Panichello
polarAngs = (category_data.Panichello.centers)'; %Polar Angles
[a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*52);
stimCols = [a;b];
stimCols_sRGB = LuvTosRGB([repelem(60, length(polarAngs)); stimCols]);
colvals = im2double(stimCols_sRGB);

p1_err_neg = category_data.Panichello.centers - category_data.Panichello.confidence(1:2:end);
p1_err_pos = mod(category_data.Panichello.confidence(2:2:end) - category_data.Panichello.centers, 360);
for i = 1:length(p1_err_neg)
errorbar(category_data.Panichello.centers(i)- offset, repelem(2,1),p1_err_neg(i), p1_err_pos(i), 'horizontal', 'ko','MarkerEdgeColor',colvals(i,:),'MarkerFaceColor',colvals(i,:));
end

% Bae et al.
polarAngs = (category_data.Bae.centers)'; %Polar Angles
% white_point = xyYToXYZ([0.3118, 0.3119, 48.64]');
% white_point = white_point ./ white_point(2);
% bae_xyz = LuvToXYZ(polarAngs,white_point);
% bae_lab = XYZtoLab(bae_xyz, white_point);
% bae_degree = rad2deg(cart2pol(bae_cielab(2,:),bae_cielab(3,:)));
%[a,b] = pol2cart(deg2rad(bae_degree),ones(1,length(bae_degree))*38);
[a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*38);
stimCols = [a;b];
stimCols_sRGB = LuvTosRGB([repelem(70, length(polarAngs)); stimCols]);
colvals = im2double(stimCols_sRGB);

p0_err_neg = category_data.Bae.centers - category_data.Bae.confidence(1:2:end);
p0_err_pos = mod(category_data.Bae.confidence(2:2:end) - category_data.Bae.centers, 360);
for i = 1:length(p0_err_neg)
errorbar(category_data.Bae.centers(i) - offset, repelem(1,1),p0_err_neg(i), p0_err_pos(i),'horizontal', 'ko','MarkerEdgeColor',colvals(i,:),'MarkerFaceColor',colvals(i,:));
end


ylim([0.5 7.5])
yticks([1:7])
% yticklabels({'Buster', 'Castor', 'Pollux', 'Humans (pooled)', 'Panichello et al.', 'Bae et al.'})
yticklabels({'Bae et al.', 'Panichello et al.', 'Humans (pooled)', 'M4', 'M3', 'M2', 'M1'})
p = gca;
p.YAxis.TickLength = [0 0];
% xlim([(offset - 360) (360-offset)])
% xticks([(offset-360):90:(360-offset)])
xlim([0 360]);
xticks([0:90:360]);
