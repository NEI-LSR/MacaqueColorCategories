% For converting files from `.mat` files to `.csv` format

datafiles = dir('*.mat');

for observer = 1:4
    load(datafiles(observer).name)

    writetable(cell2table(...
        [trialdata.dirname',...
        trialdata.cues,...
        trialdata.chosen,...
        trialdata.choices],...
        'VariableNames', {'dirname', 'cues', 'chosen', 'choices'}),...
        [datafiles(observer).name(1:end-4),'.csv'])
end

