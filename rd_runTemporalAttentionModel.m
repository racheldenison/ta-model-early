% rd_runTemporalAttentionModel.m

% wrapper for rd_temporalAttentionModel

%% setup
% nTrials = 10000;
nTrials = 1;
soas = [100:50:500 800];
% soas = [100 250 800];
endoConds = {'no-endo','endoT1','endoT2','endoT1T2'};

%% run conditions
for iEndo = 1:numel(endoConds)
    endoCond = endoConds{iEndo};
    
    for iSOA = 1:numel(soas)
        soa = soas(iSOA);
        if nTrials > 1
            fprintf('\nSOA = %d', soa)
        end
        
        %% run trials
        d = []; e = [];
        for iTrial = 1:nTrials
            if mod(iTrial,100)==0, fprintf('.'), end
            [d(:,iTrial) e(:,:,iTrial)] = rd_temporalAttentionModel(soa, endoCond);
        end
        
        %% summary
        acc = d;
        acc(d==-1) = 0;
        
        % store acc mean
        accMean(:,iSOA,iEndo) = mean(acc,2);
        
        evMean(:,iSOA,iEndo) = mean(e(:,end,:),3);
    end
end
fprintf('\n')

%% resort endo condition data into valid, invalid, neutral
if isequal(endoConds,{'no-endo','endoT1','endoT2','endoT1T2'})
    evValidity{1}(1,:) = evMean(1,:,2); % T1 valid
    evValidity{1}(2,:) = evMean(1,:,3); % T1 invalid
    evValidity{1}(3,:) = evMean(1,:,4); % T1 neutral - endoT1T2
    evValidity{1}(4,:) = evMean(1,:,1); % T1 neutral - no endo
    
    evValidity{2}(1,:) = evMean(2,:,3); % T2 valid
    evValidity{2}(2,:) = evMean(2,:,2); % T2 invalid
    evValidity{2}(3,:) = evMean(2,:,4); % T2 neutral - endoT1T2
    evValidity{2}(4,:) = evMean(2,:,1); % T2 neutral - no endo
    
    % cuing effect and average across cue validities
    for iT = 1:numel(evValidity)
        evCueEff{iT} = squeeze(evValidity{iT}(1,:,:) - evValidity{iT}(2,:,:));
        evCueAve{iT} = squeeze(mean(evValidity{iT},1));
    end
end

%% plot figs
if numel(size(evMean))==2 % if 2-dimensional
    figure
    plot(repmat(soas,2,1)', accMean','.-')
    xlabel('SOA (ms)')
    ylabel('accuracy')
    ylim([0 1])
    
    figure
    plot(repmat(soas,2,1)', evMean','.-')
    xlabel('SOA (ms)')
    ylabel('evidence')
    ylim([0 2.8])
else
    % sorted by validity
    cueValidityNames = {'valid','invalid','both','none'};
    intervalNames = {'early','late'};
    evLims = [0 4];
    soaLims = [soas(1)-100 soas(end)+100];
    colors = get(0,'DefaultAxesColorOrder');
    axTitle = '';
    
    figure
    for iT = 1:numel(evValidity)
        subplot(1,numel(evValidity),iT)
        hold on
        %     plot(soaLims, [0.5 0.5], '--k');
        
        p1 = plot(repmat(soas',1,numel(cueValidityNames)),...
            evValidity{iT}', '.-', 'MarkerSize', 20);
        
        xlabel('soa')
        ylabel('evidence')
        title(intervalNames{iT})
        xlim(soaLims)
        ylim(evLims)
        
        if iT==1
            legend(p1, cueValidityNames,'location','south')
        end
        %     rd_supertitle(subjectID);
        %     rd_raiseAxis(gca);
        %     rd_supertitle(axTitle);
    end
    
    figure
    for iT = 1:numel(evValidity)
        subplot(1,numel(evValidity),iT)
        hold on
        
        plot(soas, evCueEff{iT})
        plot(soas, evCueAve{iT},'k')
        
        xlabel('soa')
        ylabel('cuing effect / average evidence')
        title(intervalNames{iT})
        xlim(soaLims)
        ylim(evLims)
    end
    legend('cuing effect','average evidence')
end
