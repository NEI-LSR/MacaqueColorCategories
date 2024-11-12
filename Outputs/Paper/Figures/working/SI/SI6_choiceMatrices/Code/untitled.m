clear,clc, close all

load matlab.mat

try
    choice_counts_diag = choice_counts'; % transposing to match similarity matrix (so cue on x-axis, choice on y-axis)
catch
    % choiceProb_diag = model{1,1}.choice_probability;
    error('model variable structure is nested') % TODO work out why
end

for i = 1:size(choice_counts_diag,1)
    choice_counts_diag(i,:) = circshift(choice_counts_diag(i,:),i-(size(choice_counts_diag,1))/2);
end

plotSimilarityMatrix(choice_counts_diag,'choice_counts','../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)


%%

try
    presentation_counts_diag = presentation_counts'; % transposing to match similarity matrix (so cue on x-axis, choice on y-axis)
catch
    % choiceProb_diag = model{1,1}.choice_probability;
    error('model variable structure is nested') % TODO work out why
end

for i = 1:size(presentation_counts_diag,1)
    presentation_counts_diag(i,:) = circshift(presentation_counts_diag(i,:),i-(size(presentation_counts_diag,1))/2);
end

plotSimilarityMatrix(presentation_counts_diag,'presentation_counts','../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)

% plotSimilarityMatrix(presentation_counts,'presentation_counts','../',[],false) % using the same function, but note that this is *not* a similarity matrix (that would take into account the specific interactions between the available choices on each trial)
