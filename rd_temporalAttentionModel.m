function [decision evidence] = rd_temporalAttentionModel(soa, endoCond)
%
% Notes on settings for specific models:
% 1) Original + extended
%   IOR is on
%   extendResponse = 1;
%   evidenceCeiling = 1;
%   integratorType = '1-stage'
%   ceiling = 1.6
%
% 2) Limited within a span
%   distributeEndo = 1;
%   evidenceCeiling = 1;
%   integratorType = '1-stage'
%   ceiling = 2.3
%
% 3) Late-stage competition (aka "normalize then integrate extended     responses)
%   evidenceCeiling = 1;
%   integratorType = '2-stage'
%   ceiling = 375
%
% Note evidenceCeiling can be turned off if looking at accuracy (with noise).

%% Setup
plotFigs = 0;

% soa
if ~exist('soa','var')
    soa = 800;
end

% endo condition
if ~exist('endoCond','var')
    endoCond = 'endoT1T2'; % options: 'no-endo','endoT1','endoT2','endoT1T2'
end

% endo type
endoType = 'Gaussian'; % options: 'Gaussian', 'SSHat'

% Ex kernel shape
ExShape = 'square'; % options: 'Gaussian', 'square'

% proporational allocation of attention (based on cue validity)
% to make all or none, set prop = [1 0]
prop = [.6 .4]; % [high low], or [valid invalid]
neutralPropOp = 'max'; % 'max','mean'

% attentional normalization
normalizeAttentionFields = 0;

% compressive nonlinearity
compressiveNonlinearity = 0;

% extend response (early)
extendResponse = 0;

% evidence ceiling
evidenceCeiling = 1;

% distribute endo attention (limited allocation within a span)
distributeEndo = 0;
if distributeEndo
    span = 500; % ms
    totalAttn = 1 + soa/span;
    if totalAttn>2
        totalAttn = 2;
    end
    if strcmp(endoCond, 'endoT1T2')
        prop = distributeAttention(totalAttn, 0);
    else
        prop = distributeAttention(totalAttn, 1, prop(1));
    end
end

% integration settings (model-specific)
integratorType = '2-stage';
ceiling = 375; % 2.3 (limited) % 1.6 (original+extended) % 375 (late competition)

% evidence accumulation
noiseSigma = 0; % 0 % 4
evidenceScale = 0.02; % 0.02
bounds = [1.8 -1.8];
% bounds = [275 -275];
decisionType = 'firstCrossing'; % 'firstCrossing','endPoint'

%% Stimuli and attention components
% basic
x = 0:2000; % 0:2000
ExWidth = 10;
baselineMod = 0;
sigma = 1e-6;
contrast = 1;

% stim
% soa = 250;
stimCenters = 100 + soa*(0:1); % start at 100 so we have all positive times
nStim = numel(stimCenters);
stimWidth = 50; % 5, 50
stimAmps = ones(1,nStim).*contrast;
stimShape = 'square'; % 'gaussian','square'
stimulus = rd_nmMakeStim(x, stimCenters, stimWidth, stimAmps, stimShape);

% exo attention
Ax = stimCenters + 100 - round(stimWidth/2); % 100
AxWidth = 30;
Apeak =  2;
Abase = 1;
attnGain = NaN;
attnGainX = rd_nmMakeStim(x, Ax, AxWidth, repmat(1,1,nStim), 'gaussian');
% attnGainX = NaN; Ax = NaN;

% endo attention
if strcmp(endoCond,'no-endo')
    EndoGain = NaN;
else
    if strcmp(endoType,'Gaussian')
        endoShift = 30; % can shift endo to align with the response instead of the stimulus
    else
        endoShift = 0;
    end
    Endox = stimCenters + endoShift;
    switch endoCond
        case 'endoT1'
            endoProps = [prop(1) prop(2)]; % 1 = high, 2 = low
        case 'endoT2'
            endoProps = [prop(2) prop(1)];
        case 'endoT1T2'
            switch neutralPropOp
                case 'mean'
                    endoProps = [mean(prop) mean(prop)];
                case 'max'
                    endoProps = [max(prop) max(prop)];
            end
        otherwise
            error('endoCond not recognized')
    end
    EndoxWidth = AxWidth*2; % *2 *200/30
    EndoAmps = repmat(Apeak-Abase,1,numel(Endox)).*endoProps;
    EndoGain = rd_nmMakeStim(x, Endox, EndoxWidth, EndoAmps, 'gaussian');
    
    switch endoType
        case 'SSHat'
            % surround suppression (Mexican-hat style)
            SSxShift = 80;
            SSx = Endox + SSxShift;
            SSxWidthCenter = AxWidth*3;
            SSAmpsCenter = EndoAmps;
            SSGainCenter = rd_nmMakeStim(x, SSx, SSxWidthCenter, SSAmpsCenter, 'gaussian');
            SSGainSurround = rd_nmMakeStim(x, SSx, SSxWidthCenter*2, SSAmpsCenter/2, 'gaussian');
            SSGain = SSGainCenter - SSGainSurround;
            % flatten top
            SSGainCutoff = 0.7*(max(SSGain)-min(SSGain)) + min(SSGain);
            SSGain(SSGain>SSGainCutoff) = SSGainCutoff;
            EndoGain = SSGain;
    end
end

% IOR
% IORx = Ax + 200;
IORx = stimCenters + 300 - round(stimWidth/2);
IORxWidth = AxWidth*4;
IORAmps = repmat((Apeak-Abase)/2, 1, nStim);
% IORGain = rd_nmMakeStim(x, IORx, IORxWidth, IORAmps, 'gaussian');
IORGain = NaN;

% symmetrical suppression
ISx = stimCenters;
ISxWidth = AxWidth*4;
ISAmps = repmat(Apeak-Abase,1,numel(ISx));
% ISGain = rd_nmMakeStim(x, ISx, ISxWidth, ISAmps, 'gaussian');
ISGain = NaN;

% Alternatively, derive inhibitory attention component from excitatory one
IAGain = NaN; % set to nan to turn off convolution
IAxWidth = IORxWidth;
IAxAmp = IORAmps(1);
IAxKernel = makeGaussian(x,0,IAxWidth);
% crop and shift kernel
IAxShift = 0;
IAxKernel = IAxKernel(IAxKernel>max(IAxKernel)/50);
IAxKernel = [IAxKernel(end-IAxShift+1:end) IAxKernel(1:end-IAxShift)];

%% Attention field
if isnan(Ax)
    attnGain = ones(size(stimulus));
    attnGainX = ones(size(x));
elseif ~isnan(attnGain)
    % do nothing - just use attnGain as given
elseif ~isnan(attnGainX)
    attnGain = (Apeak-Abase)*attnGainX + Abase;
else
    attnGainX = makeGaussian(x,Ax,AxWidth,1);
    attnGain = (Apeak-Abase)*attnGainX + Abase;
end

%% Add Endo
if ~isnan(EndoGain)
    attnGain = attnGain + EndoGain; %%%% +
end

%% Add IOR or other inhibitory component
if ~isnan(IORGain)
    attnGain = attnGain - IORGain; %%%% -
elseif ~isnan(ISGain)
    attnGain = attnGain - ISGain; %%%% -
elseif ~isnan(IAGain)
    % convolve with excitatory attention gain
    IA = conv(attnGain,IAxKernel);
    IA = IA(1:length(attnGain));
    attnGain = attnGain - IA + IAxAmp;
end

%% Normalize attention fields
if normalizeAttentionFields
    AIxWidth = AxWidth*6; % *6
    AIxKernel = makeGaussian(x(1:1000),500,AIxWidth); % note, the kernel has to be no longer than it needs to be
    % Suppressive drive
    AI = conv2sepYcirc(attnGain - Abase, AIxKernel);
    
    % Normalization
    attnGain0 = attnGain; % just store original attnGain (for debugging)
    attnGain = attnGain0 ./ (AI + 1 + sigma);
    
    if plotFigs
        figure
%         figure(1)
        hold on
        plot(x,attnGain0)
        plot(x,AI,'k')
        plot(x,attnGain,'g')
        legend('orig attn gain','attentional suppressive drive','resulting attn gain')
    end
end

%% Stimulation field and suppressive field
switch ExShape
    case 'Gaussian'
        % ExKernel = makeGaussian(x,0,ExWidth);
        ExKernel = makeGaussian(x,100,ExWidth);
    case 'square'
        ExKernel = makeSquareWave(x,100,ExWidth);
end

% Crop kernels
ExKernel = ExKernel(ExKernel>max(ExKernel)/50);

%% Stimulus drive
% Eraw = conv2sepYcirc(stimulus,ExKernel) + baselineMod;
% We don't want circular convolution for time domain
Eraw = conv(stimulus,ExKernel) + baselineMod;
Eraw = Eraw(1:size(stimulus,2));
Emax = max(Eraw(:));
E = attnGain .* Eraw;

%% Compressive nonlinearity
if compressiveNonlinearity
    logistic = @(x,L,k,x0) L./(1 + exp(-k.*(x-x0)));
    
    E0 = E;
    E = logistic(E0,1,25,1); % 1,25,1
%     % test the logistic (for debugging)
%     x = 0:.01:2;
%     y = logistic(x,1,15,1);
%     plot(x,y)

    if plotFigs
        figure
%         figure(2)
        hold on
        plot(x,E0)
        plot(x,E,'r')
        legend('orig E','resulting E')
    end
end

%% Extend response
if extendResponse
    expk = exp(-.1*(0:.1:20)); % exponential kernel % -.02
    c = conv(E,expk);
    c1 = c*sum(E)/sum(c);
    E = c1(1:length(x));
end

%% Accumulate evidence for the decision
decisionWindowDur = min(diff(stimCenters));
switch integratorType
    case '1-stage'
        noise = noiseSigma*randn(size(x));
        % evidence = cumsum(E) + noise; % additive decision noise
        % evidence = cumsum(E + noise); % additive sensory noise
        
        % going to need separate decisions for the two stim
        % assume that the accumulation starts at stimulus onset and lasts until the
        % next stimulus. the second stimulus will have the same accumulation time
        % as the first. or can set it separately.
        stimStarts = round(stimCenters - stimWidth/2);
        for iStim = 1:nStim
            decisionWindows(iStim,:) = stimStarts(iStim):stimStarts(iStim)+decisionWindowDur;
        end
        
        evidence = zeros(size(decisionWindows));
        for iStim = 1:nStim
            dw = decisionWindows(iStim,:);
            evidence(iStim,:) = cumsum(E(dw) + noise(dw)); % additive sensory noise
        end
        evidence = evidence*evidenceScale;
        
        % evidence ceiling
        if evidenceCeiling
            evidence(evidence>ceiling) = ceiling;
        end
        
    case '2-stage'
        sensoryNoise = repmat(noiseSigma*randn(size(x)),nStim,1); % same sensory noise applied to all stim
%         sensoryNoise = noiseSigma*randn(nStim,length(x)); % different sensory noise for different stim
%         decisionNoise = noiseSigma/100*randn(nStim,length(x));
        
        stimStarts = round(stimCenters - stimWidth/2);
        dwSelector = zeros(nStim, length(x));
        
        for iStim = 1:nStim
            decisionWindows(iStim,:) = stimStarts(iStim):stimStarts(iStim)+decisionWindowDur;
            dwSelector(iStim, decisionWindows(iStim,:)) = 1;
            % evidence for each stim should listen only to that stim
            Estim(iStim,:) = E.*dwSelector(iStim,:);
            
            spreadResponse = 1;
            if spreadResponse
                expk = exp(-.02*(0:.1:60)); % exponential kernel % -.01
                Estim0 = Estim;
                c = conv(Estim(iStim,:),expk);
                c1 = c*sum(Estim(iStim,:))/sum(c);
                Estim(iStim,:) = c1(1:length(x));
            end
        end
        
        % normalize responses
        n = 2;
        sigma = 0.01;
        Estim0 = Estim;
        for iStim = 1:nStim
            Estim(iStim,:) = Estim0(iStim,:).^n./(sigma + sum(Estim0.^n,1));
        end
        
        % integrate
        for iStim = 1:nStim
            evidence(iStim,:) = cumsum(Estim(iStim,:) + sensoryNoise(iStim,:)); % additive sensory noise
        end
        
        % evidence ceiling
        if evidenceCeiling
            evidence(evidence>ceiling) = ceiling;
        end
        
%         evidence = integrateWithRateNormalization(x,Estim,4);
        
        % plot normalize and integrate
        if plotFigs
            figure
            subplot(3,1,1)
            plot(x,Estim0)
            title('Estim before normalization')
            subplot(3,1,2)
            plot(x,Estim)
            title('Estim after normalization')
            subplot(3,1,3)
            plot(x,evidence)
            title('integrated evidence')
        end
        
%         % normalize evidence pools
%         evidence0 = evidence + 0; % arbitrary constant added to prevent zero crossings, which do weird things when normalized
%         for iStim = 1:nStim
%             evidenceNorm(iStim,:) = evidence0(iStim,:)./(50 + sum(evidence0,1));
%         end
% %         evidence1 = evidence0.*evidenceNorm;
%         evidence1 = evidenceNorm;
% %         evidence = evidence + decisionNoise;
% 
%         % see if any of the traces has reached the bound
%         % if so, cut off the original at that point and re-run
%         % note! this is not set up right now to deal with the negative
%         % bound
%         bounds = [0.4 -0.4]; % bounds getting reset here!
%         cc = find(evidence1 > bounds(1), 1, 'first');
%         [ccStim ccTime] = ind2sub(size(evidence1),cc);
%         evidence0b = evidence0;
%         evidence0b(ccStim, ccTime:end) = evidence0b(ccStim, ccTime);
%         for iStim = 1:nStim
%             evidenceNormb(iStim,:) = evidence0b(iStim,:)./(50 + sum(evidence0b,1));
%         end
% %         evidenceb = evidence0b.*evidenceNormb;
%         
%         % set evidence
%         evidence = evidenceNormb;


%         % max out integrated evidence at bound
%         bounds = [60 -60]; % bound getting reset here!
%         evidence0 = evidence;
%         evidence(evidence>bounds(1)) = bounds(1);
%         evidence(evidence<bounds(2)) = bounds(2);
%         
%         % another method: online normalization
%         onInput = Estim;
% %         onInput = evidence;
%         R = integrateAndNormalize(x, onInput);
% %         R = onlineNormalization(x, onInput);
%         if plotFigs
%             figure
%             subplot(2,1,1)
%             plot(x,onInput)
% %             title('before online normalization')
%             title('before integrate & normalize')
%             subplot(2,1,2)
%             plot(x,R)
% %             title('after online normalization')
%             title('after integrate & normalize')
%         end
        
%         if plotFigs
%             figure
%             subplot(2,1,1)
%             plot(x,evidence0)
%             ylabel('evidence')
%             legend('T1','T2')
%             title('before normalization')
%             subplot(2,1,2)
%             plot(x,evidenceNorm)
%             ylabel('evidence')
%             legend('T1','T2')
%             title('after normalization')
% %             subplot(3,1,3)
% %             plot(x,evidence1)
% %             ylabel('evidence')
% %             legend('T1','T2')
% %             title('rate normalized?')
%             
%             figure
%             subplot(2,1,1)
%             hold on
%             plot(x,evidence0b)
% %             plot(x([1 end]), [bounds(1) bounds(1)], '--k')
%             ylabel('evidence')
%             legend('T1','T2')
%             title('before normalization, with bound')
%             subplot(2,1,2)
%             hold on
%             plot(x,evidenceNormb)
%             plot(x([1 end]), [bounds(1) bounds(1)], '--k')
%             ylabel('evidence')
%             legend('T1','T2')
%             title('after normalization, with bound')
% %             subplot(3,1,3)
% %             hold on
% %             plot(x,evidenceb)
% %             plot(x([1 end]), [bounds(1) bounds(1)], '--k')
% %             ylabel('evidence')
% %             legend('T1','T2')
% %             title('rate normalized?')
%         end
    otherwise
        error('integratorType not recognized')
end

%% Determine which decision was selected
% correct decision = 1, incorrect decision = -1
switch decisionType
    case 'firstCrossing'
        % We could take the first boundary crossing
        for iStim = 1:nStim
            cc = find(evidence(iStim,:) >= bounds(1), 1, 'first');
            ic = find(evidence(iStim,:) <= bounds(2), 1, 'first');
            
            if isempty(cc) && isempty(ic) % no crossing
                d = 0;
            elseif isempty(ic) && cc > 0  % only correct crossing
                d = 1;
            elseif isempty(cc) && ic > 0 % only incorrect crossing
                d = -1;
            else
                d = (cc < ic) - (ic < cc); % both cross, find the earliest
            end
            decision(iStim,1) = d;
        end
    case 'endPoint'
        % ... or the end point
        decision = (evidence(:,end) > repmat(bounds(1),2,1)) - ...
            (evidence(:,end) < repmat(bounds(2),2,1));
end

% if there is no decision from evidence, guess randomly
guesses = (round(rand(nStim,1))-0.5)*2;
decision(decision==0) = guesses(decision==0);
    
%% Plot figs
if plotFigs
%     figure
    figure(1)
    % cla
    hold on
    plot(x, stimulus,'k')
    plot(x, attnGain,'k')
    plot(x, Eraw,'c')
    plot(x, E,'g')
    switch integratorType
        case '1-stage'
            for iStim = 1:nStim
                dw = decisionWindows(iStim,:);
                plot(x(dw), evidence(iStim,:));
            end
        case '2-stage'
            plot(x, evidence);
    end
    plot(IORGain,'r')
    % plot(x, I,'r')
    % plot(x, R,'b')
    xlabel('time')
    ylabel('activity')
    ylim([-3 3])
    
    % legend('stim','attn','Eraw','E','I','R')
    legend('stim','attn','Eraw','E','evidence')
    
    title(sprintf('soa = %d ms, %s', soa, endoCond))
end
