function stimulus = rd_nmMakeStim(x, stimCenters, stimWidth, stimAmps, stimShape)
%
% function stimulus = rd_nmMakeStim(x, stimCenters, stimWidth, stimAmps)
%
% Make stimulus vector (1D) for NMOA
% stimCenters is a vector with the position of each center
% stimWidth is a scalar (same width for all stimuli currently)
% stimAmp (optional) is the amplitude of the Gaussians, a vector with the
%   amplitude for each center
%
% Rachel Denison
% Jan 2014

% stimCenters = [100 -100];

if nargin < 5 || isempty(stimShape)
    stimShape = 'gaussian';
end
if nargin < 4 || isempty(stimAmps)
    stimAmps = [1 1];
end

for iStim = 1:numel(stimCenters)
    switch stimShape
        case 'gaussian'
            stim(iStim,:) = makeGaussian(x,stimCenters(iStim),stimWidth,stimAmps(iStim));
        case 'square'
            stim(iStim,:) = makeSquareWave(x,stimCenters(iStim),stimWidth,stimAmps(iStim));
        otherwise
            error('stimShape not recognized')
    end
end
stimulus = sum(stim,1);