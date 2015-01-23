function [R] = integrateWithRateNormalization(time, I, bound)
% function [R] = integrateWithRateNormalization(time, I, bound)
%
% online integration, where the integration rate is normalized by the
% activity in the accumulators
% see also onlineNormalize.m, integrateAndNormalize.m
%
% I is the input (i.e. from a lower area) at each  time step
%
% Rachel Denison
% Jan 2015

if nargin<1 || isempty(time)
    time = 0:0.01:2;
end
if nargin<2 || isempty(I)
    I(1,:) = ones(size(time));
    I(2,:) = ones(size(time))*2;
end
if nargin<3 || isempty(bound)
    bound = NaN;
end

% tau = 1;
B = 1;
w = 100;
sigma = 10;

R(1:size(I,1),1) = 0;

for t = 1:length(time)-1
    if t==500
        a = 1;
    end
    % calculate normalization factor
%     k = (R(:,t) + B)./(sigma + sum(w*R(:,t)));
    k = (I(:,t) + B)./(sigma + sum(w*I(:,t)));
    
    % integrate
    R(:,t+1) = R(:,t) + k.*I(:,t);
    
    % limit to bound
    grBound = R(:,t+1)>bound;
    R(grBound,t+1) = bound;
end


