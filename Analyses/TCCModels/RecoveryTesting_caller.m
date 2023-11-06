clear, clc, close all

rns = 0:9; % random number seeds

for i = 1%:length(rns)
    [NLL(i),AIC(i),BIC(i)] = RecoveryTesting('real',rns(i));
end

%%

writetable(struct2table(NLL),...
    ['NLL','_',datestr(now,'yymmdd-HHMMSS'),'.csv'])

writetable(struct2table(AIC),...
    ['AIC','_',datestr(now,'yymmdd-HHMMSS'),'.csv'])

writetable(struct2table(BIC),...
    ['BIC','_',datestr(now,'yymmdd-HHMMSS'),'.csv'])

%%

% writematrix(NLL,...
%     ['NLL','_',datestr(now,'yymmdd-HHMMSS'),'.csv'])
% 
% writematrix(AIC,...
%     ['AIC','_',datestr(now,'yymmdd-HHMMSS'),'.csv'])
% 
% writematrix(BIC,...
%     ['BIC','_',datestr(now,'yymmdd-HHMMSS'),'.csv'])
