function [unique_difficulties,pct_correct,...
    difficulties_completed_counts, f]...
    = difficulty_psychometric(cleandata,filename)
%% PURPOSE:
%Error analysis by difficulty
%Get difficulty for each trial, plot errors by difficulty
%Difficulty calculated as distance of closest distractor to cue

%% Process Inputs
cues = cell2mat(cleandata.trialdata.cues);
choices = cell2mat(cleandata.trialdata.choices);
chosen = cell2mat(cleandata.trialdata.chosen);
paradigm = unique(cleandata.trialdata.paradigm);

% If real data, get dirname from cleandata
if isfield(cleandata.trialdata, 'dirname') == 1
    dirname = unique(cleandata.trialdata.dirname);
end

nBig = size(cleandata.trialdata.stimCols{1,1},1);
nSmall = sum(~isnan(cleandata.trialdata.choices{end,1}));


%% Calculate difficulty

completed_idx = ~isnan(chosen);
correct_idx = cues(completed_idx) == chosen(completed_idx);

distance = abs(cues(completed_idx) - choices(completed_idx,:));
distance(distance > nBig/2) = nBig - distance(distance > nBig/2);

for trial = 1:size(cues(completed_idx),1)
    try
        difficulty(trial) = min(nonzeros(distance(trial,:)));
    catch
        if ~any(distance(trial,:) ~= 0) % if there are no choices where distance ~= 0
            difficulty(trial) = 0;
        else
            error('unknown error'); % I can't predict any situation where this would be reached, but just in case...
        end
    end
end

unique_difficulties = 1:max(unique(difficulty));

for i = 1:max(unique_difficulties)
    difficulties_completed_counts(i,1) = sum(difficulty==i);
end

for i = 1:max(unique_difficulties)
    difficulties_correct_counts(i,1) = sum(difficulty(correct_idx)==i);
end

% Percent correct by distance of closest distractor
pct_correct = difficulties_correct_counts./difficulties_completed_counts;
pct_correct(isnan(pct_correct)) = 0;


%% Fit Weibull function with 4 params

Weibull = fittype(@(slope, inflect, floor, ceil, x) (floor+(1-floor-ceil).*(1-exp(-(x./slope).^inflect)))); %Define version of Weibull with 4 args

%% Weighted fit
[f, gof] = fit(unique_difficulties',pct_correct,...
    Weibull,...
    'weights', difficulties_completed_counts,...
    'Lower', [0 -20 0 0],...
    'Upper', [20 20 1 1],...
    'StartPoint', [10 0.01 0.5 0.5]);

%% FIGURES

figure(1) % Uncomment this to get them all on the same graph
% figure,

% axes1 = axes;
hold on
plot(unique_difficulties, pct_correct, '.','MarkerEdgeColor','none'); % this is invisible, yet neccessary, for reasons I don't understand
% fitplot = plot(f); set(fitplot,'color','k'); set(fitplot,'LineWidth',2);

xlim([min(unique_difficulties), max(unique_difficulties)]);
xticks([min(unique_difficulties) max(unique_difficulties)]);

plot(f,'k')

p22 = predint(f,unique_difficulties,0.95,'functional','on');
x_plot =[unique_difficulties, fliplr(unique_difficulties)]; % h/t https://www.mathworks.com/matlabcentral/answers/425206-plot-of-confidence-interval-with-fill
y_plot=[p22(:,1)', flipud(p22(:,2))'];
fill(x_plot, y_plot, 1,'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.2);

ylim([0.25 1]);
yticks([0.25 1]);
yticklabels({'25%', '100%'});
title('Accuracy by Trial Difficulty');
xlabel('Distance of Closest Distractor');
ylabel('Accuracy');
legend('off')

%% 

saveas(gcf,['../','difficulty_', filename, ...
    datestr(now,'yymmdd-HHMMSS'), '.svg'])

end
