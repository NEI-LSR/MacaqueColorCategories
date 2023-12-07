%% Recreating Panichello Figures 

%% Figure 1B
% Angular Error vs. Density by Load and Delay Time
%% Format Data
data = table2array(FIG2Bdata);

%% Human 
human_1short(:,3) = data(~isnan(data(:,1)));
human_1short(:,1) = 1;
human_1short(:,2) = 1;

idx = ~isnan(data(:,2));
human_1long(:,3) = data(idx,2);
human_1long(:,1) = 1;
human_1long(:,2) = 7;

idx = ~isnan(data(:,3));
human_3short(:,3) = data(idx,3);
human_3short(:,1) = 3;
human_3short(:,2) = 1;

idx = ~isnan(data(:,4));
human_3long(:,3) = data(idx,4);
human_3long(:,1) = 3;
human_3long(:,2) = 7;
%%
human = {human_1short; human_1long; human_3short; human_3long};

%% Monkey W
idx = ~isnan(data(:,5));
monkeyw_1short(:,3) = data(idx,5);
monkeyw_1short(:,1) = 1;
monkeyw_1short(:,2) = 0.5;

idx = ~isnan(data(:,6));
monkeyw_1long(:,3) = data(idx,6);
monkeyw_1long(:,1) = 1;
monkeyw_1long(:,2) = 1.5;

idx = ~isnan(data(:,7));
monkeyw_2short(:,3) = data(idx,7);
monkeyw_2short(:,1) = 2;
monkeyw_2short(:,2) = 0.5;

idx = ~isnan(data(:,8));
monkeyw_2long(:,3) = data(idx,8);
monkeyw_2long(:,1) = 2;
monkeyw_2long(:,2) = 1.5;
%%
monkeyw = {monkeyw_1short; monkeyw_1long; monkeyw_2short; monkeyw_2long};


%% Monkey E
idx = ~isnan(data(:,9));
monkeye_1short(:,3) = data(idx,9);
monkeye_1short(:,1) = 1;
monkeye_1short(:,2) = 0.5;

idx = ~isnan(data(:,10));
monkeye_1long(:,3) = data(idx,10);
monkeye_1long(:,1) = 1;
monkeye_1long(:,2) = 1.5;

idx = ~isnan(data(:,11));
monkeye_2short(:,3) = data(idx, 11);
monkeye_2short(:,1) = 2;
monkeye_2short(:,2) = 0.5;

idx = ~isnan(data(:,12));
monkeye_2long(:,3) = data(idx, 12);
monkeye_2long(:,1) = 2;
monkeye_2long(:,2) = 1.5;
%%
monkeye = {monkeye_1short; monkeye_1long; monkeye_2short; monkeye_2long};

%% Figure 1B
% Error vs. Density by condition

data = human;

data_1short = data{1,1};
data_1long = data{2,1};
data_3short = data{3,1};
data_3long = data{4,1};


val_1short = unique(data_1short(:,3));
val_1long = unique(data_1long(:,3));
val_3short = unique(data_3short(:,3));
val_3long = unique(data_3long(:,3));

for j = 1:length(val_1short)
    i = val_1short(j,:);
    error_count = data_1short(data_1short(:,3)==i, :); % sort all incorrect trials for cue i
    error_1short(j,:) = histc(error_count(:,3),i);
end

error_1short = error_1short ./sum(error_1short);

for j = 1:length(val_1long)
    i = val_1long(j,:);
    error_count = data_1long(data_1long(:,3)==i, :); % sort all incorrect trials for cue i
    error_1long(j,:) = histc(error_count(:,3),i);
end

error_1long = error_1long ./sum(error_1long);

for j = 1:length(val_3short)
    i = val_3short(j,:);
    error_count = data_3short(data_3short(:,3)==i, :); % sort all incorrect trials for cue i
    error_3short(j,:) = histc(error_count(:,3),i);
end

error_3short = error_3short ./sum(error_3short);

for j = 1:length(val_3long)
    i = val_3long(j,:);
    error_count = data_3long(data_3long(:,3)==i, :); % sort all incorrect trials for cue i
    error_3long(j,:) = histc(error_count(:,3),i);
end

error_3long = error_3long ./sum(error_3long);

%% Plots

% Human load 1 short delay
figure;
hold on
%bar(val_1short, error_1short);
gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d'; %a = amp, b = mu, c = stdv, d = vert. offset
%startingPoints = [0 50 22 0];
[f,gof]= fit(val_1short, error_1short, gaussEqn);
short1 = plot(f,val_1short, error_1short);
set(short1(2),'LineWidth',2);
set(short1(1),'visible','off');
[f,gof]= fit(val_1long, error_1long, gaussEqn);
long1 = plot(f,val_1long, error_1long);
set(long1(2),'Linestyle','--','LineWidth',2);
set(long1(1),'visible','off');
% [f,gof]= fit(val_3short, error_3short, gaussEqn, 'start',startingPoints);
short3 = plot(f,val_3short, error_3short);
set(short3(2),'Color','b','LineWidth',2);
set(short3(1),'visible','off');
[f,gof]= fit(val_3long, error_3long, gaussEqn);
long3 = plot(f,val_3long, error_3long);
set(long3(2),'Linestyle','--','Color','b','LineWidth',2);
set(long3(1),'visible','off');
legend('hide');
xlabel('Angular error');
ylabel('Density');
title('Monkey W');
xlim([-180 180]);
hold off

%% Calculate Error
for i = 1:length(data) % for each incorrect trial
    if abs(data(i,1) - data(i,2)) < 180
        data(i,3) = (data(i,1) - data(i,2));
    elseif abs(data(i,1) - data(i,2)) > 180 && data(i,1) > data(i,2)
    data(i,3) = (-(180 - abs(data(i,1) - data(i,2))));
    else
    data(i,3) = (180 - abs(data(i,1) - data(i,2)));  
    end
end

data = array2table(data);
data.Properties.VariableNames = {'Cue', 'Chosen', 'Error'};
%%
%monkeyw = table2array(monkeyw);
data = monkeyw
cues = data(:,1);
chosen = data(:,2);
error = data(:,3);
all_cues = unique(cues);
%%
cue_dist = [cues, chosen];
%%
for i = 1:length(all_cues)
    cue_value = all_cues(i,1);
    dist_count = cue_dist(cue_dist(:,1)==cue_value, 2); % sort all incorrect trials for cue i
    error_mean(i,:) = mean(dist_count);
    %error_counts(:,i) = histc(dist_count(:,2),total_distances);
end
 
%% Plot biases (Figure 4d)
figure('WindowState','maximized');
% human 4d
subplot(1,3,1)
hold on
plot(1:length(human4d),human4d(:,1));
plot(1:length(human4d), human4d(:,2));
hold off
legend('Empirical', 'Model');
ylim([-20 20]);
xlim([0 90]);
ylabel('Bias (mean response - target)');
xlabel('Target color');
title('Human');

% monkey W 4d
subplot(1,3,2)
hold on
plot(1:length(monkeyw4d),monkeyw4d(:,1));
plot(1:length(monkeyw4d), monkeyw4d(:,2));
hold off
legend('Empirical', 'Model');
ylim([-50 50]);
ylabel('Bias (mean response - target)');
xlabel('Target color');
title('Monkey W');

subplot(1,3,3)
hold on
plot(1:length(monkeye4d),monkeye4d(:,1));
plot(1:length(monkeye4d), monkeye4d(:,2));
hold off
legend('Empirical', 'Model');
ylim([-50 50]);
ylabel('Bias (mean response - target)');
xlabel('Target color');
title('Monkey E');
