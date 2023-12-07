function model = fitMixtureModel(cleandata,lengthOfSlidingWindow,includeCorrect)

try nBig = cleandata.nBig;
catch
    try
        nBig = size(cleandata.trialdata.allchoices{1,1},2);
    catch
        nBig = 64;
        warning('Assuming nBig = 64')
    end
end
try nSmall = cleandata.nSmall;
catch
    try
        nSmall = size(cleandata.trialdata.choices{1,1},2);
    catch
        nSmall = 4;
        warning('Assuming nSmall = 4')
    end
end
interval = 360/nBig;

%% Filter data
% Remove aborted trials, and correct trials (we WANT the incorrect trials)

if istable(cleandata)
    cues    = cleandata.cues;
    choices = [cleandata.choices_1,cleandata.choices_2,cleandata.choices_3,cleandata.choices_4];
    chosen  = cleandata.chosen;
else
    try
        cues    = cell2mat(cleandata.trialdata.cues);
        choices = cell2mat(cleandata.trialdata.choices);
        chosen  = cell2mat(cleandata.trialdata.chosen);
    catch % for Panichello data
        cues    = cleandata.trialdata.cues;
        choices = cell2mat(cleandata.trialdata.choices);
        chosen  = cleandata.trialdata.chosen;
    end
end

abortIndex = or(isnan(chosen),any(isnan(choices'))');
correctIndex = cues == chosen;
if ~exist('includeCorrect','var') || ~includeCorrect
    filter = or(abortIndex,correctIndex);
else
    filter = abortIndex;
end

cues_filtered       = cues(~filter);
choices_filtered    = choices(~filter,:);
chosen_filtered     = chosen(~filter);

nTrials_filtered = sum(~filter);

% % _No filter mode_
% cues_filtered       = cues;
% choices_filtered    = choices;
% chosen_filtered     = chosen;
%
% nTrials_filtered = length(cues);

%% Calculate Angular Error

% Distance values
PotentialDistances = (-180+interval:interval:180)';
hue_angle_by_index = 0:interval:360-interval;

% Calculate angular error (distance) between incorrect choice and cue

d = zeros(nTrials_filtered,1);
for trial = 1:nTrials_filtered
    d(trial) = rad2deg(angdiff(...
        deg2rad(chosen_filtered(trial)*interval),...
        deg2rad(cues_filtered(trial)  *interval)));
end
d_index = round(d/interval);
d_index(d_index == -nBig/2) = nBig/2;
d = d_index*interval;

% Count of number of times each color (by distance from cue) was chosen
choice_counts = zeros(length(PotentialDistances),nBig);
bin_edges = (-180+(interval/2):interval:180+(interval/2));

for i = 1:nBig
    choice_counts(:,i) = histcounts(d(cues_filtered == i), bin_edges);
end

% figure,
% imagesc(choice_counts')
% colorbar
% axis square
%
% figure,
% choice_counts_midExtract = [...
%     choice_counts(1:31,:); ...
%     NaN(1,nBig); ...
%     choice_counts(33:end,:)...
%     ]';
% imagesc(choice_counts_midExtract,...
%     'AlphaData', ~isnan(choice_counts_midExtract))
% colorbar                                                                    %TODO: Set NaN to be highlighted as a different thing on the colorbar. I'm sure I've done that before...
% axis square

%% Count of number of times each color (by distance from cue) was an option
% (completed trials only)

presentation_counts = zeros(nBig);

cues_completed    = cues(~abortIndex);
choices_completed = choices(~abortIndex,:);

for cueIndex = 1:nBig

    choices_filtered_forThisCueIndex = choices_completed(cues_completed == cueIndex,:); % choices for (completed) trials matching this cueIndex
    for choice = 1:nSmall
        for trial = 1:size(choices_filtered_forThisCueIndex,1)
            choices_filtered_forThisCueIndex(trial,choice) = rad2deg(angdiff(...
                deg2rad(hue_angle_by_index(choices_filtered_forThisCueIndex(trial,choice))),...
                deg2rad(hue_angle_by_index(cueIndex))));
            choices_filtered_forThisCueIndex(trial,choice) = round(choices_filtered_forThisCueIndex(trial,choice)/interval);
            if choices_filtered_forThisCueIndex(trial,choice) == -nBig/2
                choices_filtered_forThisCueIndex(trial,choice) = nBig/2;
            end
            choices_filtered_forThisCueIndex(trial,choice) = choices_filtered_forThisCueIndex(trial,choice)*interval;
        end
    end
    for PotentialDistanceIndex = 1:length(PotentialDistances)
        presentation_counts(PotentialDistanceIndex,cueIndex) =...
            sum(abs(choices_filtered_forThisCueIndex(:) - PotentialDistances(PotentialDistanceIndex)) < 1);
    end
end

% figure,
% imagesc(presentation_counts')
% colorbar
% axis square
%
% figure,
% presentation_counts_midExtract = [...
%     presentation_counts(1:31,:); ...
%     NaN(1,nBig); ...
%     presentation_counts(33:end,:)...
%     ]';
% imagesc(presentation_counts_midExtract,...
%     'AlphaData', ~isnan(presentation_counts_midExtract))
% colorbar                                                                    %TODO: Set NaN to be highlighted as a different thing on the colorbar. I'm sure I've done that before...
% axis square

%% Normalize counts for each color by number of times color was an option

choice_probability = choice_counts./presentation_counts;

% figure,
% imagesc(choice_probability')
% colorbar
% axis square

% Replace value at 0 (correct choice) to include from curve fit
if ~exist('includeCorrect','var') || ~includeCorrect
    choice_probability(nBig/2,:) = NaN;
    presentation_counts(nBig/2,:) = NaN;
else
end

%% Fit Gaussian to error distribution for each cue to get bias value

gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d';

startingPoints = [0.5 0 60 0.1];

for cueIndex = 1:nBig

    dist_idx = find(~isnan(choice_probability(:,cueIndex))); % Index of distance counts for only colors shown

    f = fit(PotentialDistances(dist_idx,:), choice_probability(dist_idx,cueIndex), ...
        gaussEqn, 'Weights', presentation_counts(dist_idx,cueIndex),...
        'start',startingPoints, 'Lower',[0 -180 0 0],'Upper',[Inf 180 Inf 1]);

    model.gaussfits{cueIndex} = f; % model output                                        % TODO Ideally this would happen at the model fitting stage rather than being tagged on here
    bias(cueIndex) = f.b;

    ci = confint(f,0.95);
    ci_lower_95(cueIndex,1) = ci(1,2);
    ci_upper_95(cueIndex,1) = ci(2,2);
    if any(isnan(ci),'all')
        disp(cueIndex)
        disp(f)
        disp(ci)
        warning('NaN in CI')
    end
end

%% Category center locations

if ~exist('lengthOfSlidingWindow','var') || isempty(lengthOfSlidingWindow)
    lengthOfSlidingWindow = 3; % moving average input
end

if mod(lengthOfSlidingWindow,2) == 0
    error('lengthOfSlidingWindow should be odd')
end

moving_bias = movmean(padarray(bias',lengthOfSlidingWindow,"circular"),lengthOfSlidingWindow,'Endpoints','discard');
moving_bias = moving_bias(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

be_w = moving_bias([1:end,1]); % bias estimates including wraparound

for i = 1:nBig
    crossesZero(i) = and(be_w(i)>=0, be_w(i+1)<=0);
end

crossings = find(crossesZero);

hue_angle_w = 0:interval:360; %includes wraparound

% Category Center Location Interpolation
if ~isempty(crossings)
    for i = 1:length(crossings)
        x = [hue_angle_w(crossings(i)) hue_angle_w(crossings(i)+1)];
        y = [be_w(crossings(i)) be_w(crossings(i)+1)];
        interp_crossing(i,1) = interp1(y,x,0);
    end
else
    interp_crossing = [];
end

%% Confidence intervals of category centers
if ~isempty(ci)

    lower_95 = movmean(padarray(ci_lower_95,lengthOfSlidingWindow,"circular"),...
        lengthOfSlidingWindow,'Endpoints','discard');
    lower_95 = lower_95(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

    upper_95 = movmean(padarray(ci_upper_95,lengthOfSlidingWindow,"circular"),...
        lengthOfSlidingWindow,'Endpoints','discard');
    upper_95 = upper_95(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

    for i = 1:nBig % check whether the confidence interval for each cue includes 0
        withinCI(i) = and(lower_95(i,1)<=0, upper_95(i,1)>=0);
    end

    with_CI = withinCI([1:end 1]);


    for i = 1:nBig
        changes(i) = or(and(with_CI(i)==0, with_CI(i+1)==1), and(with_CI(i)==1, with_CI(i+1)==0));
    end

    change_range = find(changes);

    upper_95_w = upper_95([1:end 1]);
    lower_95_w = lower_95([1:end 1]);

    if ~isempty(change_range) % if confidence interval spans zero for all values

        change_range = change_range([1:end 1]);

        crossings = interp_crossing / interval;

        for i = 1:length(crossings)
            for j = 1:sum(changes)
                if change_range(j) < change_range(j+1)
                    if (crossings(i) >= change_range(j) && crossings(i) <= change_range(j+1)) == 1
                        CI_range(i,:) = [change_range(j) change_range(j+1)];
                    end
                elseif change_range(j) > change_range(j+1)
                    if (crossings(i) >= change_range(j) && crossings(i) >= change_range(j+1)) == 1 ...
                            || (crossings(i) <= change_range(j) && crossings(i) <= change_range(j+1) == 1)
                        CI_range(i,:) = [change_range(j) change_range(j+1)];
                    end
                end
            end
        end

        CI_range = CI_range';
        CI_range = CI_range(:);

        % Confidence Interval Interpolation
        for i = 1:length(CI_range)
            if be_w(CI_range(i)) > 0
                x = [hue_angle_w(CI_range(i)) hue_angle_w(CI_range(i)+1)];
                y = [lower_95_w(CI_range(i)) lower_95_w(CI_range(i)+1)];
            elseif  be_w(CI_range(i)) < 0
                x = [hue_angle_w(CI_range(i)) hue_angle_w(CI_range(i)+1)];
                y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
            end

            interp_ci(i,1) = interp1(y,x,0);

            if isnan(interp_ci(i,1)) && be_w(CI_range(i)) > 0
                x = [hue_angle_w(CI_range(i)) hue_angle_w(CI_range(i)+1)];
                y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
                interp_ci(i,1) = interp1(y,x,0);
            elseif isnan(interp_ci(i,1)) && be_w(CI_range(i)) < 0
                x = [hue_angle_w(CI_range(i)) hue_angle_w(CI_range(i)+1)];
                y = [lower_95_w(CI_range(i)) lower_95_w(CI_range(i)+1)];
                interp_ci(i,1) = interp1(y,x,0);
            end
        end
        %     else
        %         %interp_ci = []; %needed?
    end
else
    interp_ci = [];
end

%% Pack up the model

model.PotentialDistances = PotentialDistances;
model.interp_crossing = interp_crossing;
model.interp_ci = interp_ci;
model.presentation_counts = presentation_counts;
model.choice_probability = choice_probability;
model.lower_95_w = lower_95_w;
model.upper_95_w = upper_95_w;
model.bias = bias;
model.be_w = be_w;
model.ci = ci;
model.moving_bias = moving_bias;

end