function [P,h] = quantPlot(X,Y,Q)
% [P,h] = quantPlot(X,Y,Q)
% draws quantiles of vector, fills them, plots median
% X is the x-axis values to be averaged over 
% Y are the values to be averaged for the Y axis
% Q are the quantiles to plot, can be scalar (e.g. 20 will plot quantiles
% 20% from the median, or a vector of quantiles (e.g. [0 15 50 85 100]).
% The default is to use quantiles corresponding to standard deviations of a
% normal distribution.
% 
% John Grogan, 2019

if ~exist('Q','var') || isempty(Q)% default percentiles - normal SD quantiles
    percs = [0 0.2 2.3 15.9 50 84.1 97.7 99.8 100];
elseif all(size(Q)==1) % if single value, use as step size in quantiles from median
        percs = [fliplr(50:-Q:0) 50:Q:100];
        percs = [0 percs 100]; % make sure has 0 and 100
        percs = unique(percs,'stable');%remove duplicates
else
    percs = Q;
end

if all(size(Y)>1)
    
    n = length(percs); % number of percentiles
    
    
    for i = 1:n % get percentiles
        P(i,:) = prctile(Y',percs(i));
    end
    
    m = median(1:n); % get median row
    
    alpha = 1; % transparency
    
    colours = {'b','g','y','r','r','y'}; % colour order
    
    for i = 1:(n-1)/2
        
        xx = [X' fliplr(X')]; % get x twice in row
        
        yy = [P(i,:) fliplr(P(n-i+1,:))]; % y coords for plotting
        
        %     remove NaNs
        xx(isnan(yy)) = [];
        yy(isnan(yy)) = [];
        
        %plot
        h(i,2) = fill(xx ,yy, colours{i}, 'linestyle','none','FaceAlpha',alpha);
        
        hold on;
        
    end
else
    P = Y';
    m=1;
end
% plot mean line on top of area
h(1,1) = plot( X', P(m,:), 'k','LineWidth',1.5);
end
