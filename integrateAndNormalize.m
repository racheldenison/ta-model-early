function [R] = integrateAndNormalize(time, I, bound)
% function [R] = integrateAndNormalize(time, I, bound)
%
% online integration + normalization
% see also onlineNormalize.m
% Louie et al. 2014 J Neurosci
%
% I is the input (i.e. from a lower area) at each  time step
%
% Rachel Denison
% Jan 2015

if nargin==0
    time = 0:0.01:2;
    I(1,:) = ones(size(time));
    I(2,:) = ones(size(time))*2;
end
if ~exist('bound','var')
    bound = NaN;
end

% tau = 1;
B = 0;
w = 0.1;
sigma = .5;

R(1:size(I,1),1) = 0;

for t = 1:length(time)-1
    % integrate
    r = R(:,t) + I(:,t);
    
    % normalize
    R(:,t+1) = (r + B)./(sigma + sum(w*r));
    
    % limit to bound
    grBound = R(:,t+1)>bound;
    R(grBound,t+1) = R(grBound,t);
end


