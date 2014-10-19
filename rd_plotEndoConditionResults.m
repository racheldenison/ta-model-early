% rd_plotEndoConditionResults.m

%% load data
endo0 = load('data/no_endo_workspace.mat');
endoT1 = load('data/endoT1_workspace.mat');
endoT2 = load('data/endoT2_workspace.mat');
endoT1T2 = load('data/endoT1T2_workspace.mat');

soas = endo0.soas;

%% organize data
for iEL = 1:2
    tAcc{iEL}(:,1) = endo0.accMean(iEL,:);
    tAcc{iEL}(:,2) = endoT1.accMean(iEL,:);
    tAcc{iEL}(:,3) = endoT2.accMean(iEL,:);
    tAcc{iEL}(:,4) = endoT1T2.accMean(iEL,:);
end

% calculate cuing effect (valid - invalid)
tCueEff(:,1) = tAcc{1}(:,2) - tAcc{1}(:,3); % T1
tCueEff(:,2) = tAcc{2}(:,3) - tAcc{2}(:,2); % T2

%% plot figs
% accuracy for T1 and T2
ylims = [.75 1];
figure
for iEL = 1:2
    subplot(1,2,iEL)
    plot(soas, tAcc{iEL})
    xlabel('soa')
    ylabel('acc')
    ylim(ylims)
end
legend('no endo','endo T1','endo T2','endo T1 & T2')

% cuing effect
figure
plot(soas, tCueEff)
ylim([0 .2])
xlabel('soa')
ylabel('cueing effect (valid - invalid)')
legend('T1','T2')


