function MemToolbox2DFallonPosteriors
% MemToolbox2DFallonPosteriors
% plot heatmaps of pars against each other for two conds for 1D and 2D

% clear
% close all

% mean params from Fallon2018
k = [17.2 19 25.2 25];
kStd = [1 1.4 1.5 1.7]; % actually SEM
a = [.828 .89 .942 .94];
aStd = [.014 .005 .004 .003];
g = [.084 .072 .021 .035];
gStd = [.01 .008 .003 .005];
b = [.088 .035 .037 .025];
bStd = [.008 .005 .006 .003];

params = [g;b;k]';
params(:,3) = k2deg(params(:,3)); % convert to degrees

%%
nSamples = 50000;

numTrials = 64;
itemsPerTrial = ones(4, numTrials) .* [4;2;4;2];


for i = 1:size(params,1)
    % 1D
    data1 = GenerateDisplays(numTrials, itemsPerTrial(i,:));
    data1.errors = SampleFromModel(SwapModel(), params(i,:),[1 numTrials],data1);
    mcmc(i,1) = MCMC(data1, SwapModel(),'PostConvergenceSamples',nSamples);
    % PlotPosterior(mcmc(1), {'g','b','SD'});

    % 2D
    data2 = GenerateDisplays2D(numTrials, itemsPerTrial(i,:));
    data2.errors = SampleFromModel2D(SwapModel2D(), params(i,:),[1 numTrials],data2);

    mcmc(i,2) = MCMC(data2, SwapModel2D(),'PostConvergenceSamples',nSamples);
    % PlotPosterior(mcmc(2), {'g','b','SD'});
end

mcmc2 = mcmc;
%% transform data

% b and g use logit/arcsine/log
mcmc = mcmc2;
for j = 1:2
    for i = 1:size(params,1)
%         mcmc(i,j).vals(:,1:3) = log(mcmc2(i,j).vals(:,1:3)); % log
%         mcmc(i,j).vals(:,1:2) = log( mcmc2(i,j).vals(:,1:2)./(1-mcmc2(i,j).vals(:,1:2)) ); % logit
%         mcmc(i,j).vals(:,1:2) = norminv(mcmc2(i,j).vals(:,1:2)); % probit
        mcmc(i,j).vals(:,1:2) = asin( sqrt( mcmc2(i,j).vals(:,1:2) ) ); % arcsine

    end
end
%%
modelNames = {'1D';'2D'};
f = figure(1);clf;
inds = [1 2; 1 3; 2 3;];
parNames = {'\gamma', '\beta','SD'};
cmaps = colormap('hsv');
colours = ['r','y','g','b'];
condInds = [1 2 3 4];
condNames = {'Ig','T1','Up','T2'};
n=13;
for j = 1:2
    for i = 1:3
        subplot(2,3,(j-1)*3+i)
        hold on
        
        for k = 1:length(condInds)
            [V,C] = hist3(mcmc(condInds(k),j).vals(:,[inds(i,:)]), [15 15]);
            V = V ./ max(V(:));
            ch=imagesc(C{1}, C{2}, V' + (n*(k-1)));
            set(ch, 'AlphaData', V','CDataMapping','direct');
            h = 1/9*ones(3);
            VS = filter2(h,V);
            [~,maxInd] = max(mcmc(condInds(k),j).like);
%             map(k,:,j) = mcmc(condInds(k),j).vals(maxInd,inds(i,:));
%             p(k,i,j) = plot(map(k,1,j), map(k,2,j),'.',...
%                 'Color',colours(k),'MarkerSize',10);
            p(k,i,j) = plot(NaN, NaN,'s','Color',colours(k),...
                'MarkerFaceColor',colours(k),'MarkerSize',4); % plot NaN so colours show up in legend
%             plot(params(k,inds(i,1)), params(k,inds(i,2)),'^',...
%                 'Color',colours(k),'MarkerSize',6,'MarkerFaceColor',colours(k));
        end
        
        box off
        axis tight;
        box off
        xlabel(parNames{inds(i,1)},'FontWeight','bold')
        if i==1
            ylabel(sprintf('%s\n\n%s',modelNames{j},parNames{inds(i,2)}),'FontWeight','bold')
        else
            ylabel(parNames{inds(i,2)},'FontWeight','bold')
        end

    end
end
legend(condNames(condInds),'Location','SouthEast')
% 
makeSubplotScalesEqual(2,3,[1 4])
makeSubplotScalesEqual(2,3,[2 5])
makeSubplotScalesEqual(2,3,[3 6])
% colorbar('Ticks',1:n:n*length(condInds),'TickLabels',condNames(condInds),'Limits',[1 n*length(condInds)],'Location','east');

% saveas(figure(1),'./Figs/MemToolbox2DFallonPosteriors.jpg')
%%

save('./Data/MemToolbox2DFallonPosteriors.mat','mcmc','params','nSamples','numTrials','itemsPerTrial','mcmc2');

end