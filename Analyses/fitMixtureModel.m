function [moving_bias, lower_95, upper_95] = fitMixtureModel(cleandata,Lab,lengthOfSlidingWindow)

if isfield(cleandata.trialdata,'dirname') % for real data
    nBig = size(cleandata.trialdata.stimCols{1,1},1);
else
    nBig = size(cleandata.trialdata.stimCols,2); % for sim data
end

nSmall = sum(~isnan(cleandata.trialdata.choices{end,1}));
interval = 360/nBig;

%% Filter data
% Remove aborted trials, and correct trials (we WANT the incorrect trials)

cues    = cell2mat(cleandata.trialdata.cues);
choices = cell2mat(cleandata.trialdata.choices);
chosen  = cell2mat(cleandata.trialdata.chosen);

abortIndex = isnan(chosen);                                                 % Note: this previously looked at cues, choices, and chosen, for NaNs, not just at choices. If you have problems, this might be why.
correctIndex = cues == chosen;
filter = or(abortIndex,correctIndex);

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
PotentialDistances = (-180+interval:interval:180)';                         % TODO: Double check the logic that this married with "32" being treated as the zero point in the rest of this function

% Calculate angular error (distance) between incorrect choice and cue

d = zeros(nTrials_filtered,1); 
for trial = 1:nTrials_filtered
    d(trial) = rad2deg(angdiff(deg2rad(chosen_filtered(trial)*interval), deg2rad(cues_filtered(trial)*interval)));
end

d = round(d,4); % cut off number of digits after decimal point

% Count of number of times each color (by distance from cue) was chosen
choice_counts = zeros(length(PotentialDistances),nBig);
bin_edges = (-180+interval:interval:180+interval);

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

% for cueIndex = 1:nBig
%     comp_trial = choices_filtered(cues_filtered == cueIndex,:); % choices for (completed) trials matching this cueIndex
%     for choice = 1:nSmall
%         for trial = 1:size(comp_trial,1)
%             if      abs(comp_trial(trial,choice) - cueIndex) < nBig/2           
%                 comp_trial(trial,choice) = (comp_trial(trial,choice) - cueIndex) * interval;
%             elseif  abs(comp_trial(trial,choice) - cueIndex) > nBig/2 && comp_trial(trial,choice) > cueIndex
%                 comp_trial(trial,choice) = (-(nBig - abs(comp_trial(trial,choice) - cueIndex))) * interval;
%             else
%                 comp_trial(trial,choice) = (nBig - abs(comp_trial(trial,choice) - cueIndex)) * interval;
%             end
%         end
%     end
%     for m = 1:length(PotentialDistances)
%         presentation_counts(m,cueIndex) = sum(comp_trial(:) == PotentialDistances(m));
%     end
% end

for cueIndex = 1:nBig
    comp_trial = choices_filtered(cues_filtered == cueIndex,:); % choices for (completed) trials matching this cueIndex
    for choice = 1:nSmall
        for trial = 1:size(comp_trial,1)
            comp_trial(trial,choice) = rad2deg(angdiff(deg2rad(comp_trial(trial,choice)*interval), deg2rad(cueIndex*interval)));
            comp_trial = round(comp_trial,4);
        end
    end
    for m = 1:length(PotentialDistances)
       presentation_counts(m,cueIndex) = sum(comp_trial(:) == PotentialDistances(m));
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

% Replace value at 0 (correct choice) to exclude from curve fit
choice_probability(nBig/2,:) = NaN;
presentation_counts(nBig/2,:) = NaN; 

%% Fit Gaussian to error distribution for each cue to get bias value

gaussEqn = 'a*exp(-(((x-b)^2)/(2*c^2)))+d';

startingPoints = [0.5 0 60 0.1];

for cueIndex = 1:nBig  

    dist_idx = find(~isnan(choice_probability(:,cueIndex))); % Index of distance counts for only colors shown

    f = fit(PotentialDistances(dist_idx,:), choice_probability(dist_idx,cueIndex), ...
        gaussEqn, 'Weights', presentation_counts(dist_idx,cueIndex),...
        'start',startingPoints, 'Lower',[0 -180 0 0],'Upper',[Inf 180 Inf 1]);
    
    mo{cueIndex} = f; % model output                                        % TODO Ideally this would happen at the model fitting stage rather than being tagged on here

    % f = fit(PotentialDistances, choice_probability(:,cueIndex), ...                       % ...I tried, and failed, here
    %     gaussEqn, 'Weights', presentation_counts(:, cueIndex),...
    %     'start',startingPoints, 'Lower',[0 -180 0 0],'Upper',[Inf 180 Inf 1]);
    % 

    bias(cueIndex) = f.b;

    ci = confint(f,0.95);
    % ci_lower_95(cueIndex) = ci(1,2);
    % ci_upper_95(cueIndex) = ci(2,2);
    if any(isnan(ci),'all')
        disp(cueIndex)
        disp(f)
        disp(ci)
        warning('NaN in CI')
    end

    % temporary plot (TODO will move into plotting later)
    figure('visible','off')
    hold on
    axis tight
    ylim([0,1])

    pltCol = 1 - repmat(...
        (presentation_counts(:,i)-min(presentation_counts(:,i)))...
        /(max(presentation_counts(:,i))-min(presentation_counts(:,i))),...
        1,3);

    s = scatter(PotentialDistances, choice_probability(:,cueIndex),'filled');
    s.CData = pltCol;
    
    % plot(f,PotentialDistances,choice_probability(:,cueIndex),'k.')
    plot(f,dists,choice_probability(~isnan(choice_probability(:,cueIndex)),cueIndex),'k.');
    
    p = gca;
    p.Children(2).Marker = 'none'; % turn off data, so that we can replot it how we like...
    p.Children(1).LineWidth = 3;

    p11 = predint(f,dists,0.95,'functional','off');                         %TODO: Check whether this is the appropriate type of interval: https://www.mathworks.com/help/curvefit/confidence-and-prediction-bounds.html
    plot(dists,p11,'k:',...
        'DisplayName','Nonsimultaneous Functional Bounds')
    p.Children(1).LineWidth = 1;
    p.Children(2).LineWidth = 1;

    xline(0,'k--')
    
    xlabel('Error distance')
    ylabel('Choice Probability')
    text(-90,0.9,num2str(cueIndex))
    legend('off')
    drawnow

    saveas(gcf,fullfile([num2str(cueIndex),'_mixMod_BreakOut.svg']))
    close all
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
    [a,s] = pol2cart(deg2rad(polarAngs),ones(1,length(polarAngs))*37);
    crossingCols = [a;s];
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


% Stuff to return...
%
% - dirname