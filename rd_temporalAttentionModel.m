function [decision evidence] = rd_temporalAttentionModel(soa, endoCond)

%% Setup
plotFigs = 1;

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

% attentional normalization
normalizeAttentionFields = 1;

%% Stimuli and attention components
% basic
x = 0:2000; % 0:2000
ExWidth = 10;
baselineMod = 0;
sigma = 1e-6;

% stim
% soa = 250;
stimCenters = 100 + soa*(0:1); % start at 100 so we have all positive times
nStim = numel(stimCenters);
stimWidth = 50; % 5, 50
stimAmps = ones(1,nStim);
stimShape = 'square'; % 'gaussian','square'
stimulus = rd_nmMakeStim(x, stimCenters, stimWidth, stimAmps, stimShape);

% exo attention
Ax = stimCenters + 100 - round(stimWidth/2);
AxWidth = 30;
Apeak =  2;
Abase = 1;
attnGain = NaN;
% attnGainX = rd_nmMakeStim(x, Ax, AxWidth, repmat(1,1,nStim), 'gaussian');
attnGainX = NaN; Ax = NaN;

% endo attention
if strcmp(endoCond,'no-endo')
    EndoGain = NaN;
else
    switch endoCond
        case 'endoT1'
            Endox = stimCenters(1);
        case 'endoT2'
            Endox = stimCenters(2);
        case 'endoT1T2'
            Endox = stimCenters;
        otherwise
            error('endoCond not recognized')
    end
    EndoxWidth = AxWidth*2; % *2 *200/30
    EndoAmps = repmat(Apeak-Abase,1,numel(Endox));
    EndoGain = rd_nmMakeStim(x, Endox, EndoxWidth, EndoAmps, 'gaussian');
    
    switch endoType
        case 'SSHat'
            % surround suppression (Mexican-hat style)
            SSxShift = 80;
            SSx = Endox + SSxShift;
            SSxWidthCenter = AxWidth*3;
            SSAmpsCenter = repmat(Apeak-Abase,1,numel(SSx));
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

% evidence accumulation
noiseSigma = 0; % 4
evidenceScale = 0.02; % 0.02
bounds = [.8 -.8];
decisionWindowDur = min(diff(stimCenters));
decisionType = 'firstCrossing'; % 'firstCrossing','endPoint'

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
    AIxWidth = AxWidth*6;
    AIxKernel = makeGaussian(x(1:1000),500,AIxWidth); % note, the kernel has to be no longer than it needs to be
    % Suppressive drive
    AI = conv2sepYcirc(attnGain - Abase, AIxKernel);
%     AI = zeros(size(x));
%     for i = 1:numel(Endox)
%         AI = AI + rd_nmMakeStim(x, Endox(i), AIxWidth, EndoAmps(i)/2, 'gaussian') + 1;
%     end
    
    % Normalization
    attnGain0 = attnGain; % just store original attnGain (for debugging)
    attnGain = attnGain0 ./ (AI + 1 + sigma);
    
    if plotFigs
        figure(1)
        hold on
        plot(x,attnGain0)
        plot(x,AI,'k')
        plot(x,attnGain,'g')
        legend('orig attn gain','attentional suppressive drive','resulting attn gain')
    end
end

%% Stimulation field and suppressive field
% ExKernel = makeGaussian(x,0,ExWidth);
ExKernel = makeGaussian(x,100,ExWidth);

% Crop kernels
ExKernel = ExKernel(ExKernel>max(ExKernel)/50);

%% Stimulus drive
% Eraw = conv2sepYcirc(stimulus,ExKernel) + baselineMod;
% We don't want circular convolution for time domain
Eraw = conv(stimulus,ExKernel) + baselineMod;
Eraw = Eraw(1:size(stimulus,2));
Emax = max(Eraw(:));
E = attnGain .* Eraw;

%% Accumulate evidence for the decision
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

% allEv(:,:,iRun) = evidence;
% iRun = iRun + 1;

%% Determine which decision was selected
% correct decision = 1, incorrect decision = -1
switch decisionType
    case 'firstCrossing'
        % We could take the first boundary crossing
        for iStim = 1:nStim
            cc = find(evidence(iStim,:) > bounds(1), 1, 'first');
            ic = find(evidence(iStim,:) < bounds(2), 1, 'first');
            
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
    figure(2)
    % cla
    hold on
    plot(x, stimulus,'k')
    plot(x, attnGain,'k')
    plot(x, Eraw,'c')
    plot(x, E,'g')
    for iStim = 1:nStim
        dw = decisionWindows(iStim,:);
        plot(x(dw), evidence(iStim,:));
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
