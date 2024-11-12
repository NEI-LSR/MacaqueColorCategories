% input(['Warning: you are about to re-fit all the models.',newline,...
%     'This will take roughly 45 minutes.',newline,...
%     'Press enter to continue.'])

% realOrSimData = 'real';
% 
% % rn = 0;
% for rn = 1:9
%     RecoveryTesting(realOrSimData,rn)
% end

%%

plotbar_NLL_AIC_BIC('Combined','../')

%% Castor

plotbar_NLL_AIC_BIC('Castor','../')