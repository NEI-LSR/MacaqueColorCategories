function [moving_bias, lower_95, upper_95] = fitMixtureModel(cleandata,Lab,lengthOfSlidingWindow)


if isfield(cleandata.trialdata,'dirname') % for real data
    dirname = unique(cleandata.trialdata.dirname);
    paradigm = unique(cleandata.trialdata.paradigm);
    if length(paradigm) > 1
        warning('More than one paradigm in data file');
    end
    try
        nBig = size(cleandata.trialdata.allchoices{1,1},2);
    catch
        nBig = size(cleandata.trialdata.stimCols{1,1},1);
    end
else
    dirname = 'simdata'; % for simdata
    nBig = size(cleandata.trialdata.stimCols,2);
end

nSmall = sum(~isnan(cleandata.trialdata.choices{end,1}));

interval = 360/nBig;

% Colors
stimCols = generateStimCols('nBig',nBig,'sat',37);
rotVal = 360/(nBig*2);
rotationMatrix = [cosd(rotVal),-sind(rotVal);sind(rotVal),cosd(rotVal)];
stimCols_biased = rotationMatrix * stimCols;
if Lab %CIELAB
    stimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LabTosRGB([repelem(76.0693, nBig); stimCols_biased]);
else % CIELUV
    stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols]);
    rstimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols_biased]);
end
colvals = im2double(stimCols_sRGB);

%% Calculate bias 

cues = cell2mat(cleandata.trialdata.cues);
choices = cell2mat(cleandata.trialdata.choices);
chosen = cell2mat(cleandata.trialdata.chosen);

try
    colresp = [cues, chosen, choices];
catch
    colresp = [cues, chosen', choices];
end

% Get completed trials + incorrect trials
colresp(any(isnan(colresp),2),:) = []; % gets rid of all aborts
incorrect = colresp(colresp(:,1) ~= colresp(:,2),1:2); % cue + response on incorrect trials


% Calculate Angular Error

% Distance values
distances = (-180+interval:interval:180)';

% Calculate angular error (distance) between incorrect choice and cue
distance = zeros(size(incorrect(:,1)));
for i = 1:length(incorrect) % for each trial
    if abs(incorrect(i,2) - incorrect(i,1)) < nBig/2
        distance(i,1) = (incorrect(i,2) - incorrect(i,1)) * interval;
    elseif abs(incorrect(i,2) - incorrect(i,1)) > nBig/2 && incorrect(i,2) > incorrect(i,1)
        distance(i,1) = (-(nBig - abs(incorrect(i,2) - incorrect(i,1)))) * interval;
    else
        distance(i,1) = (nBig - abs(incorrect(i,2) - incorrect(i,1))) * interval;
    end
end

% Cue + angular error for each incorrect trial
cue_dist = [incorrect(:,1), distance];

% Count of number of times each color (by distance from cue) was chosen
distance_counts = zeros(length(distances),nBig);
for i = 1:nBig
    dist_count = cue_dist(cue_dist(:,1)==i, :); % All angular error values for cue i
    distance_counts(:,i) = histc(dist_count(:,2),distances);
end

% Count of number of times each color (by distance from cue) was an option
% (completed trials only)
for i = 1:nBig
    comp_trial = colresp(colresp(:,1)==i,3:end); % choices for completed trials
    for j = 1:nSmall %
        for k = 1:size(comp_trial,1)
            if abs(comp_trial(k,j) - i) < nBig/2
                comp_trial(k,j) = (comp_trial(k,j) - i)*interval;
            elseif abs(comp_trial(k,j) - i) > nBig/2 && comp_trial(k,j) > i
                comp_trial(k,j) = (-(nBig - abs(comp_trial(k,j) - i)))*interval;
            else
                comp_trial(k,j) = (nBig - abs(comp_trial(k,j) - i))*interval;
            end
        end
    end
    for m = 1:length(distances)
        compchoice_distances(m,i)= sum(comp_trial(:)==distances(m,1));
    end
end

% Normalize counts for each color by number of times color was an option
dist_prop = distance_counts./compchoice_distances;

% Replace value at 0 (correct choice) to exclude from curve fit
dist_prop(nBig/2,:) = NaN;
compchoice_distances(nBig/2,:) = NaN; % Replace number at distance = 0

%Fit Gaussian to error distribution for each cue to get bias value
for i = 1:nBig
    gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d';
    startingPoints = [0.5 0 78 0];
    dist_idx = find(~isnan(dist_prop(:,i))); % Index of distance counts for only colors shown
    dists = distances(dist_idx,:); % Excludes distance values for colors never shown
    weights = compchoice_distances(dist_idx,i); % Weights fit by number of times each color was an option
    f = fit(dists,dist_prop(~isnan(dist_prop(:,i)),i),gaussEqn,'start',startingPoints, 'Weights',weights,'Lower',[0 -180 0 0],'Upper',[Inf 180 Inf 1]);
    ci = confint(f,0.95);
    bias(i,1) = f.b;
    ci_lower_95(i,1) = ci(1,2);
    ci_upper_95(i,1) = ci(2,2);
    if any(isnan(ci),'all')
        disp(i)
        disp(f)
        disp(ci)
        warning('NaN in CI')
    end
end

%% Category center locations

if ~exist('lengthOfSlidingWindow','var')
    lengthOfSlidingWindow = 3; % moving average input
end

if mod(lengthOfSlidingWindow,2) == 0
    error('lengthOfSlidingWindow should be odd')
end

moving_bias = movmean(padarray(bias,lengthOfSlidingWindow,"circular"),lengthOfSlidingWindow,'Endpoints','discard');
moving_bias = moving_bias(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

be_w = moving_bias([1:end,1]); % bias estimates including wraparound

for i = 1:nBig
    crossesZero(i) = and(be_w(i)>=0, be_w(i+1)<=0);
end

crossings = find(crossesZero);

hue_angle = 0:interval:360; %includes wraparound

% Category Center Location Interpolation
if ~isempty(crossings)

    for i = 1:length(crossings)
        x = [hue_angle(crossings(i)) hue_angle(crossings(i)+1)];
        y = [be_w(crossings(i)) be_w(crossings(i)+1)];
        interp_crossing(i,1) = interp1(y,x,0);
    end

    % Category Center Colors
    polarAngs = interp_crossing'; %Polar Angles
    [a,b] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
    crossingCols = [a;b];
    if Lab % CIELAB
        crossingCols_sRGB = LabTosRGB([repelem(76.0693, length(polarAngs)); crossingCols]);
    else % CIELUV
        crossingCols_sRGB = LuvTosRGB([repelem(76.0693, length(polarAngs)); crossingCols]);
    end
    crossing_colvals = im2double(crossingCols_sRGB);
else
    interp_crossing = [];
end

%% Confidence intervals of category centers
if isempty(ci) == 0

    lower_95 = movmean(padarray(ci_lower_95,lengthOfSlidingWindow,"circular"),lengthOfSlidingWindow,'Endpoints','discard');
    lower_95 = lower_95(ceil(lengthOfSlidingWindow/2)+1:end-ceil(lengthOfSlidingWindow/2));

    upper_95 = movmean(padarray(ci_upper_95,lengthOfSlidingWindow,"circular"),lengthOfSlidingWindow,'Endpoints','discard');
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

    if isempty(change_range) == 0 % if confidence interval spans zero for all values

        change_range = change_range([1:end 1]);

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
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
                y = [lower_95_w(CI_range(i)) lower_95_w(CI_range(i)+1)];
            elseif  be_w(CI_range(i)) < 0
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
                y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
            end

            interp_ci(i,1) = interp1(y,x,0);

            if isnan(interp_ci(i,1)) && be_w(CI_range(i)) > 0
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
                y = [upper_95_w(CI_range(i)) upper_95_w(CI_range(i)+1)];
                interp_ci(i,1) = interp1(y,x,0);
            elseif isnan(interp_ci(i,1)) && be_w(CI_range(i)) < 0
                x = [hue_angle(CI_range(i)) hue_angle(CI_range(i)+1)];
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

end