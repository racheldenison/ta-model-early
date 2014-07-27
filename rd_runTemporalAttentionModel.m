% rd_runTemporalAttentionModel.m

% wrapper for rd_temporalAttentionModel

%% setup
nTrials = 10000;
soas = [100 150 200 250 300 400 500 700 1000];
% soas = 100;

%% run conditions
for iSOA = 1:numel(soas)
    soa = soas(iSOA);
    fprintf('\nSOA = %d', soa)
    
    %% run trials
    for iTrial = 1:nTrials
        if mod(iTrial,100)==0, fprintf('.'), end
        d(:,iTrial) = rd_temporalAttentionModel(soa);
    end
    
    %% summary
    acc = d;
    acc(d==-1) = 0;
    
    % store acc mean
    accMean(:,iSOA) = mean(acc,2);
end

%% plot figs
figure
plot(repmat(soas,2,1)', accMean','.-')
xlabel('SOA (ms)')
ylabel('accuracy')