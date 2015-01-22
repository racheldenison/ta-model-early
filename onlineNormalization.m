function [R, G] = onlineNormalization(time, V, bound)
% function [R, G] = onlineNormalization(time, V, bound)
%
% online normalization
% Louie et al. 2014 J Neurosci
%
% Rachel Denison
% Jan 2015

if nargin==0
    time = 0:0.01:2;
    V(1,:) = ones(size(time));
    V(2,:) = ones(size(time))*2;
end
if ~exist('bound','var')
    bound = NaN;
end

% tau = 1;
B = 0;
w = 1;

R(1:size(V,1),1) = 0;
G(1:size(V,1),1) = 0;

for t = 1:length(time)-1
    % excitatory
    deltaR = -R(:,t) + (V(:,t) + B)./(1 + G(:,t));
    
    % inhibitory
    deltaG = -G(:,t) + sum(w*R(:,t));
    
    % update
    R(:,t+1) = R(:,t) + deltaR;
    G(:,t+1) = G(:,t) + deltaG;
    
    % limit to bound
    grBound = R(:,t+1)>bound;
    R(grBound,t+1) = R(grBound,t);
end


