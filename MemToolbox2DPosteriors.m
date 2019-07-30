function MemToolbox2DPosterios()
% MemToolbox2DPosteriors
% MCMC posterior samples, plotted against each other, for 1D and 2D models

% clear
% close all
%%
nSamples = 10000;
params = [.3 .3 30];
numTrials = 100;
itemsPerTrial = ones(1, numTrials) .* 3;
% 1D
data1 = GenerateDisplays(numTrials, itemsPerTrial);
data1.errors = SampleFromModel(SwapModel(), params,[1 numTrials],data1);
mcmc(1) = MCMC(data1, SwapModel(),'PostConvergenceSamples',nSamples);
% PlotPosterior(mcmc(1), {'g','b','SD'});

% 2D
data2 = GenerateDisplays2D(numTrials, itemsPerTrial);
data2.errors = SampleFromModel2D(SwapModel2D(), params,[1 numTrials],data2);

mcmc(2) = MCMC(data2, SwapModel2D(),'PostConvergenceSamples',nSamples);
% PlotPosterior(mcmc(2), {'g','b','SD'});

%%
modelNames = {'1D';'2D'};
set(0,'DefaultAxesFontWeight','bold')
f = figure();
inds = [1 2; 1 3; 2 3;];
parNames = {'\gamma', '\beta','SD'};
for j = 1:2
    for i = 1:3
        subplot(2,3,(j-1)*3+i)
        [V,C] = hist3(mcmc(j).vals(:,[inds(i,:)]), [15 15]);
        imagesc(C{1}, C{2}, V');
        set(gca,'YDir','normal');
        set(gca, 'box', 'off');
        axis tight;
        box off
        xlabel(parNames{inds(i,1)})
        if i==1
            ylabel(sprintf('%s\n\n%s',modelNames{j},parNames{inds(i,2)}))
        else
            ylabel(parNames{inds(i,2)})
        end

        [correls(j,i,1),correls(j,i,2)] = corr(mcmc(j).vals(:,inds(i,1)), mcmc(j).vals(:,inds(i,2)), 'type', 'Spearman');
    end
end
colormap(palettablecolormap('sequential'));
makepalettable(f);

disp(correls)

%%

save('./Data/MemToolbox2DPosteriors.mat', 'mcmc')
end