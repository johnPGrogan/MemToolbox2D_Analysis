function MemToolbox2DPaperFigs
% MemToolbox2DPaperFigs


%% 

load('./Data/MemToolbox2DSim.mat','fitPars','params2','cols','sd','g','b');

%% guess
legendLab = {'y=x','median','±1SD','±2SD','±3SD','±4SD'};
figure(1);clf
Q = [];
n = length(cols);
parNames = {'\gamma','\beta','\sigma'};
inds = [];
for iB = 1:n
    for iSD=1:n
        if any(~isnan(col(fitPars(:,1,cols(iB),cols(iSD),:))))

            subplot(n,n,(iB-1)*n+iSD)
            p = quantPlot(sq(params2(:, 1,cols(iB), cols(iSD),1)), sq(fitPars(:,1,cols(iB),cols(iSD),:)));
            hold on
            line([0 1],[0 1],'Color','k','LineStyle','--')
            if iB==1
                title(sprintf('\\sigma=%.1f',sd(cols(iSD))),'FontWeight','normal')
            end
            if iSD==1
                ylabel(sprintf('\\beta=%.1f',b(cols(iB))))
            end
            inds = [inds, (iB-1)*n+iSD];
        end
    end
end
h = findobj(gca, 'Type', 'Line','-or','Type','Patch');
legend(h,legendLab,'Location',[.8064 .1242 .1482 .2012]);

makeSubplotScalesEqual(n,n,inds)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('Simulated \gamma value','FontWeight','normal')
h.YLabel.Visible = 'on';
y = ylabel('Recovered \gamma value','FontWeight','normal');
y.Position(1) = y.Position(1) *1.2;

saveas(figure(1), 'Figs/MemToolbox2D_gamma.jpg');

%% misbind
figure(2);clf
inds = [];
for iG = 1:n
    for iSD=1:n
        if any(~isnan(col(fitPars(cols(iG),2,:,cols(iSD),:))))
            subplot(n,n,(iG-1)*n+iSD)
            p = quantPlot(sq(params2(cols(iG), 2,:, cols(iSD),1)), sq(fitPars(cols(iG),2,:,cols(iSD),:)));
            hold on
            line([0 1],[0 1],'Color','k','LineStyle','--')
            if iG==1
                title(sprintf('\\sigma=%.0f',sd(cols(iSD))),'FontWeight','normal')
            end
            if iSD==1
                ylabel(sprintf('\\gamma=%.1f',g(cols(iG))))
            end
            inds = [inds, (iG-1)*n+iSD];
        end
    end
end
h = findobj(gca, 'Type', 'Line','-or','Type','Patch');
legend(h,legendLab,'Location',[.8064 .1242 .1482 .2012]);
makeSubplotScalesEqual(n,n,inds)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('simulated \beta value','FontWeight','bold')
h.YLabel.Visible = 'on';
y = ylabel('recovered \beta value','FontWeight','bold');
y.Position(1) = y.Position(1) *1.2;

saveas(figure(2), 'Figs/MemToolbox2D_beta.jpg');


%% SD
figure(3);clf
inds = [];
for iG = 1:n
    for iB=1:n
        if any(~isnan(col(fitPars(cols(iG),3,cols(iB),:,:))))
            subplot(n,n,(iG-1)*n+iB)
            p = quantPlot(sq(params2(cols(iG), 3, cols(iB), :,1)), sq(fitPars(cols(iG),3,cols(iB),:,:)));
            hold on
            line([0 100],[0 100],'Color','k','LineStyle','--')
            if iG==1
                title(sprintf('\\beta=%.1f',b(cols(iB))),'FontWeight','normal')
            end
            if iB==1
                ylabel(sprintf('\\gamma=%.1f',g(cols(iG))))
            end
            axis([0 100 0 100])
            inds = [inds, (iG-1)*n+iB];
        end
    end
end
h = findobj(gca, 'Type', 'Line','-or','Type','Patch');
legend(h,legendLab,'Location',[.7475 .1099 .1482 .2012]);

makeSubplotScalesEqual(n,n,inds)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('simulated \sigma value','FontWeight','bold')
h.YLabel.Visible = 'on';
y = ylabel('recovered \sigma value','FontWeight','bold');
y.Position(1) = y.Position(1) *1.2;

saveas(figure(3), 'Figs/MemToolbox2D_sigma.jpg');


%%

%% plot diff in pars vs pars
f = load('./Data/MemToolbox1DSim.mat', 'fitParsOrdered','simParsOrdered');
load('./Data/MemToolbox2DSim.mat','fitParsOrdered','simParsOrdered')
modelNames = {'1D', '2D'};
parNames = {'\sigma','\alpha','\beta','\gamma'};
fitPars1DOrdered = f.fitParsOrdered;
simPars1DOrdered = f.simParsOrdered;
d = nancat(4, fitPars1DOrdered - simPars1DOrdered, fitParsOrdered - simParsOrdered);

figure(4);clf;
simParsOrdered1(:,:,1) = round(simPars1DOrdered,2,'significant');
simParsOrdered1(:,:,2) = round(simParsOrdered,2,'significant');
simParsOrdered1(simParsOrdered1<0) = NaN;
c  = [1 0 0; 1 .4 0; 1 1 0;0 1 0; 0 1 1; 0 0 1;];
colormap(c);
for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       diff = reshape(groupMeans(d(:,i,:,j),1, repmat(simParsOrdered1(:,i,j),1,1,size(d,4)),'dim'),11,[]);
       h = errorBarPlot(diff','area',1,'standardError',2);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==1; title(parNames{i}); end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(sprintf('%s\nrecovery error',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i,j), 0+(i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([-100 100]); set(gca, 'YTick', -80:80:80);
       else; ylim([-.6 .6]); set(gca, 'YTick', -.5:.5:.5);
       end
    end
end
makeSubplotScalesEqual(2,4,[2:4, 6:8])
makeSubplotScalesEqual(2,4,[1 5])
colormap(c)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(4), 'Figs/MemToolbox2D_recovErr.jpg');

%%
clear
load('./Data/MemToolbox2DPosteriors.mat', 'mcmc')
modelNames = {'1D';'2D'};
set(0,'DefaultAxesFontWeight','bold')
f = figure(5);
inds = [1 2; 1 3; 2 3;];
parNames = {'\gamma', '\beta','\sigma'};
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
        if i==1; ylabel(sprintf('%s\n\n%s',modelNames{j},parNames{inds(i,2)}))
        else; ylabel(parNames{inds(i,2)}); end
        if j==1; title([parNames{inds(i,1)} ' & ' parNames{inds(i,2)}]);end
        [correls(j,i,1),correls(j,i,2)] = corr(mcmc(j).vals(:,inds(i,1)), mcmc(j).vals(:,inds(i,2)), 'type', 'Spearman');
        xlim([0 .65]);
        set(gca,'XTick',0:.3:.6);
        if i>1
            ylim([0 65]);
            set(gca,'YTick',0:30:60);
        else
            ylim([0 .65]);
            set(gca,'YTick',0:.3:.6);
        end
    end
end
colormap(palettablecolormap('sequential'));
makepalettable(f);

disp(correls)
saveas(figure(5),'Figs/MemToolbox2D_posteriors.jpg')

%%
clear
load('./Data/MemToolbox2DSimNTrials.mat','simParsOrdered','fitParsOrdered',...
    'fitPars1DOrdered','numTrials','modelNames')

% plot diff in pars vs pars
d = nancat(4, fitPars1DOrdered - simParsOrdered, fitParsOrdered - simParsOrdered);
parNames = {'\sigma','\alpha','\beta','\gamma'};
figure(6);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c = colormap('jet');
c = c(1:3:end,:);
for j=1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,:,j),1,simParsOrdered1(:,i,1),'dim'),[4,3,2,1]);
       h = errorBarPlot(diff(:,:,1:20,i),'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==1;title(parNames{i});end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(sprintf('%s\nrecovery error',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 0 + (i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([-40 20]); set(gca,'YTick',-40:20:20);
       else; ylim([-.4 .2]); set(gca,'YTick', -.4:.2:.2);end
    end
end
subplot(2,4,8), colorbar(gca,'South','Ticks',[0 1], 'TickLabels',[numTrials([1 end])])

makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])

h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(6), 'Figs/MemToolbox2D_nTrials.jpg');

%% plot abs diff in pars vs pars
figure(7);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c = colormap('jet');
c = c(1:3:end,:);
for j=1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,:,j),1,simParsOrdered1(:,i,1),'dim'),[4,3,2,1]);
       h = errorBarPlot(abs(diff(:,:,1:20,i)),'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==1;title(parNames{i});end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(sprintf('%s\nabs( recovery error )',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 0+(i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([0 100]); set(gca,'YTick',0:50:100);
       else; ylim([0 .5]); set(gca,'YTick', 0:.25:.5);end
    end
end
subplot(2,4,8), colorbar(gca,'North','Ticks',[0 1], 'TickLabels',[numTrials([1 end])])



makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(7), 'Figs/MemToolbox2D_nTrialsAbs.jpg');


%%
clear
load('./Data/MemToolbox2DSimNDist.mat','simParsOrdered','fitParsOrdered',...
    'fitPars1DOrdered','modelNames','nMaxDist','parNames')
% plot diff in pars vs pars
d = nancat(4, fitPars1DOrdered - simParsOrdered, fitParsOrdered - simParsOrdered);

figure(8);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [1 0 0; 1 .4 0; 1 1 0;0 1 0; 0 1 1; 0 0 1;];
for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,:,j),1,simParsOrdered1(:,i),'dim'),[4,3,2,1]);
       h = errorBarPlot(diff(:,:,:,i),'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==2,    xlabel(parNames{i},'FontWeight','bold'),
       else; title(parNames{i}); end
       if i==1, ylabel(sprintf('%s\nrecovery error',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 0+(i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([-40 20]); set(gca,'YTick',-40:20:20);
       else; ylim([-.4 .2]); set(gca,'YTick', -.4:.2:.2);end
    end
end
makeSubplotScalesEqual(2,4,[2:4, 6:8])
makeSubplotScalesEqual(2,4,[1 5])
colormap(c)
subplot(2,4,8), colorbar(gca,'South','Ticks',[0 1], 'TickLabels',[1 nMaxDist])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(8), 'Figs/MemToolbox2D_nDist.jpg');

%% 
clear
load('./Data/MemToolbox2DSimConstBias.mat','fitParsOrdered', 'simParsOrdered', ...
    'fitParsOrdered3','fitPars2','fitPars3','bias','modelNames','parNames')
labels = {'Biased','Corrected'};

d = nancat(4, fitParsOrdered - simParsOrdered, fitParsOrdered3(:,1:4,:) - simParsOrdered);

figure(9);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [ 0 0 1; 0 1 1;0 1 0; 1 1 0; 1 .4 0;1 0 0;];
for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,:,j),1,simParsOrdered1(:,i),'dim'),[4,3,2,1]);
       h = errorBarPlot(diff(:,:,:,i),'area',1);
       hold on;
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==1; title(parNames{i}); end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(sprintf('%s\nrecovery error',biasNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 0 + (i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x([1 ceil(length(x)/2) end]))
       xlim([1 11])
       if i==1; ylim([-200 200]); set(gca,'YTick',-200:100:200);
       else; ylim([-.6 .6]); set(gca,'YTick', -.5:.5:.5);end
    end
end
makeSubplotScalesEqual(2,4,[2:4, 6:8])
makeSubplotScalesEqual(2,4,[1 5])
colormap(c)
subplot(2,4,8), colorbar(gca,'South','Ticks',[0 1], 'TickLabels',arrayfun(@num2str, bias([1 end]),'UniformOutput',0))
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(9), 'Figs/MemToolbox2D_constBias.jpg');

%%

clear
load('./Data/MemToolbox2DSimPropBias.mat','fitParsOrdered','fitPars3Ordered',...
    'simParsOrdered','nSteps','bias','nParSets','modelNames','parNames')
labels = {'Uncorrected','Corrected'};

d = nancat(5, fitParsOrdered(:,1:4,:) - simParsOrdered, fitPars3Ordered(:,1:4,:) - simParsOrdered);
parNames = {'\sigma','\alpha','\beta','\gamma'};

figure(10);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c = colormap('jet');
c = c(1:11:end,:);
colormap(c);
for j=1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,:,:,j),1,simParsOrdered1(:,i),'dim'),[4,3,2,1]);
       h = errorBarPlot(diff(:,:,:,i),'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if j==1, title(parNames{i}), end
       box off
       x = unique(round(simParsOrdered1(:,i), 1));
       set(gca,'XTick',1:4:11,'XTickLabel',x([1 ceil(length(x)/2) end]))
       xlim([1 11])
       if i==1
           ylabel(sprintf('%s\nrecovery error',biasNames{j}),'FontWeight','bold')
       end
       if i==1; ylim([-100 100]); set(gca,'YTick',linspace(-80,80,3));
       else; ylim([-.25 .25]); set(gca,'YTick', linspace(-.2,.2,3));end;
    end
end
subplot(2,4,5), 
h = colorbar(gca,'North','Ticks',[0 .5 1],'TickLabels',[bias(1) 0 bias(end)]);
h.Position = [0.3114 0.5048 0.3939 0.0279];
makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('simulated parameter','FontWeight','normal')
x.Position(2) = x.Position(2)*1.2;

saveas(figure(10), 'Figs/MemToolbox2D_propBias.jpg');



%%

clear
load('./Data/MemToolbox2DSimRadBias.mat','fitParsOrdered','fitPars3Ordered',...
    'simParsOrdered','nSteps','bias','nParSets','modelNames','parNames')
labels = {'Uncorrected','Corrected'};

d = nancat(5, fitParsOrdered(:,1:4,:) - simParsOrdered, fitPars3Ordered(:,1:4,:) - simParsOrdered);
parNames = {'\sigma','\alpha','\beta','\gamma'};


figure(11);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [ 0 0 1; 0 1 1;0 1 0; 1 1 0; 1 .4 0;1 0 0;];
colormap(c);
for j=1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,:,:,j),1,simParsOrdered1(:,i),'dim'),[4,3,2,1]);
       h = errorBarPlot(diff(:,:,:,i),'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if j==1, title(parNames{i}), end
       if i==1, ylabel(sprintf('%s\nrecovery error',biasNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 0 + (i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x([1 ceil(length(x)/2) end]))
       xlim([1 11])
       if i==1; ylim([-80 80]); set(gca,'YTick',-50:50:50);
       else; ylim([-.1 .1]); set(gca,'YTick', -.05:.05:.05);end;
    end
end
subplot(2,4,5), colorbar(gca,'North','Ticks',[0 .5 1],'TickLabels',[bias(1) 0 bias(end)])
makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(11), 'Figs/MemToolbox2D_radBias.jpg');


%%

clear

load('./Data/MemToolbox2DSimDistrRadial.mat','fitParsOrdered','fitPars3Ordered','simParsOrdered','bias','nIter')
labels = {'Uncorrected','Corrected'};
% plot diff in pars vs pars
d = nancat(5, fitParsOrdered(:,1:3,:) - simParsOrdered, fitPars3Ordered(:,1:3,:) - simParsOrdered);
parNames = {'\sigma','\alpha','\gamma'};


figure(12);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [ 0 0 1; 0 1 1;0 1 0; 1 1 0; 1 .4 0;1 0 0;];
colormap(c);
for j=1:2
    for i = 1:3
       subplot(2,3,(j-1)*3+i)
       set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,:,:,j),1,simParsOrdered1(:,i),'dim'),[4,3,2,1]);
       h = errorBarPlot(diff(:,:,:,i),'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if j==1, title(parNames{i}), end
       if i==1; ylabel(sprintf('%s\nrecovery error',biasNames{j}),'FontWeight','bold');end
       box off
       x = unique(round(simParsOrdered1(:,i), 0 + (i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x([1 ceil(length(x)/2) end]))
       xlim([1 11])
       if i==1; ylim([-120 120]); set(gca,'YTick',-100:100:100);
       else; ylim([-.75 .75]); set(gca,'YTick', -.5:.5:.5);end
    end
end
subplot(2,3,4), colorbar(gca,'North','Ticks',[0 .5 1],'TickLabels',[bias(1) 0 bias(end)])
makeSubplotScalesEqual(2,3,[2:3 5:6])
makeSubplotScalesEqual(2,3,[1 4])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2) .* 1.2;

saveas(figure(12), 'Figs/MemToolbox2D_distrRad.jpg');

%%
clear

load('./Data/MemToolbox2DSimStimConstr.mat','fitParsOrdered','simParsOrdered','parNames')
% plot diff in pars vs pars
d = cat(4, fitParsOrdered - simParsOrdered, abs(fitParsOrdered - simParsOrdered));
for i = 1:4
    [~,diffP(i)] = ttest(d(:,i,1), d(:,i,2));
end
labels = {'Recovery error','Abs( recovery error )'};
modelNames = {'Unconstrained','Constrained'};
figure(13);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [1 0 0; 1 .4 0; 1 1 0;0 1 0; 0 1 1; 0 0 1;];
for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
%        set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,i,:,j),1,simParsOrdered1(:,i,:),'dim'),[3,1,2]);
       h = errorBarPlot(diff,'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==2; xlabel(parNames{i},'FontWeight','bold');
       else; title(parNames{i}); end
       if i==1, ylabel(labels{j},'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i,1), 0+(i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([-18 18]); set(gca,'YTick',-15:15:15);
       else; ylim([-.06 .06]); set(gca,'YTick', -.05:.05:.05);end
    end
end
makeSubplotScalesEqual(2,4,[2:4])
makeSubplotScalesEqual(2,4,[6:8])
legend(h(:,1),modelNames,'Location','Best')
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;
saveas(figure(13), 'Figs/MemToolbox2D_stimConstr.jpg');

%%

clear
load('./Data/MemToolbox2DSimResampling.mat','fitParsOrdered','simParsOrdered')
parNames = {'\sigma','\alpha','\beta','\gamma'};
% plot diff in pars vs pars
d = cat(4, fitParsOrdered - simParsOrdered, abs(fitParsOrdered - simParsOrdered));
for i = 1:4
    [~,diffP(i)] = ttest(d(:,i,1), d(:,i,2));
end
labels = {'Recovery error','Abs( recovery error )'};
modelNames = {'Resampling','Constraining'};
figure(14);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [1 0 0; 1 .4 0; 1 1 0;0 1 0; 0 1 1; 0 0 1;];
for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
%        set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,i,:,j),1,simParsOrdered1(:,i,:),'dim'),[3,1,2]);
       h = errorBarPlot(diff,'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==1, title(parNames{i}); end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(labels{j},'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i,1), 0+(i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([-18 18]); set(gca,'YTick',linspace(-15, 15, 3));
       else; ylim([-.06 .06]); set(gca,'YTick', linspace(-.05, .05, 3));end
    end
end
makeSubplotScalesEqual(2,4,[2:4])
makeSubplotScalesEqual(2,4,[6:8])
legend(h(:,1),modelNames,'Location','Best')
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;
saveas(figure(14), 'Figs/MemToolbox2D_resampling.jpg');

%%
clear
load('./Data/MemToolbox2DSim2AFC.mat','fitParsOrdered','simParsOrdered','modelNames')
parNames = {'\sigma','\alpha','\gamma'};
% plot diff in pars vs pars
d = fitParsOrdered - simParsOrdered;

figure(15);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [1 0 0; 1 .4 0; 1 1 0;0 1 0; 0 1 1; 0 0 1;];
for j = 1:2
    for i = 1:3
       subplot(2,3,(j-1)*3+i)
%        set(gca,'ColorOrder',c);
       diff = permute(groupMeans(d(:,:,j),1,simParsOrdered1(:,i),'dim'),[3,2,1]);
       h = errorBarPlot(diff(:,:,i),'area',1);
       hold on
       line([0 100],[0 0],'Color','k','LineStyle','--')
       if j==1; title(parNames{i}); end
       if j==2, xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(sprintf('%s\nrecovery error',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 0 + (i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([-80 80]); set(gca,'YTick',linspace(-60,60, 3));
       else; ylim([-.3 .3]); set(gca,'YTick', linspace(-.25,.25, 3));end
    end
end
makeSubplotScalesEqual(2,3,[2:3, 5:6])
makeSubplotScalesEqual(2,3,[1 4])
colormap(c)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2) .* 1.2;

saveas(figure(15), 'Figs/MemToolbox2D_2AFC.jpg');


%%
clear
load('./Data/MemToolbox2DFallonPosteriors.mat');

modelNames = {'1D';'2D'};
f = figure(16);clf;
inds = [1 2; 1 3; 2 3;];
parNames = {'\gamma', '\beta','\sigma'};
cmaps = colormap(hsv(4));
colormap(cmaps);

condInds = [1 2 3 4];
condNames = {'Ignore','T1','Update','T2'};
n=13;
for j = 1:2
    for i = 1:3
        subplot(2,3,(j-1)*3+i)
        hold on
        
        for k = 1:length(condInds)
            [V,C] = hist3(mcmc(condInds(k),j).vals(:,[inds(i,:)]), [20 20]);
            V = V ./ max(V(:));
            ch=imagesc(C{1}, C{2}, V' + (n*(k-1)));
            set(ch, 'AlphaData', V','CDataMapping','scaled');
            h = 1/9*ones(3);
            VS = filter2(h,V);
            [~,maxInd] = max(mcmc(condInds(k),j).like);
            p(k,i,j) = plot(NaN, NaN,'s','Color',cmaps(k,:),...
                'MarkerFaceColor',cmaps(k,:),'MarkerSize',4); % plot NaN so colours show up in legend
        end
  
        box off
        axis tight;
        box off
        if j==1
            title([parNames2{inds(i,1)} ' & ' parNames2{inds(i,2)}],'FontWeight','bold');
        else
            xlabel(parNames2{inds(i,1)},'FontWeight','bold')
        end
        if i==1
            ylabel(sprintf('%s\n\n%s',modelNames{j},parNames2{inds(i,2)}),'FontWeight','bold')
            ylim([0 .6]);
            set(gca,'YTick',0:.3:.6);
        else
            ylabel(parNames2{inds(i,2)},'FontWeight','bold')
            ylim([7 25]);
            set(gca,'YTick',10:10:30);
        end
        xlim([0 .6]);
        set(gca,'XTick', 0:.3:.6);
    end
end
legend(condNames(condInds),'Location','Best')
% 
makeSubplotScalesEqual(2,3,[1 4])
makeSubplotScalesEqual(2,3,[2 5])
makeSubplotScalesEqual(2,3,[3 6])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2) .* 1.2;

saveas(figure(16),'./Figs/MemToolbox2D_fallonPosteriors.jpg')


%%

clear
load('./Data/MemToolbox2DBehMetrics.mat');
varNames = {'targDist', 'nearestNeighbour','targDist - nearestNeighbour','swapError','swapErrorCorr','swapErrorMeanThresh','swapErrorMeanThreshCorr'};
parNames = {'\sigma','\alpha','\beta','\gamma'};

figure(17);
inds = [1 2; 1 4; 2 1; 3 3; 4 2];
correls4 = NaN(2,5,2);
for i = 1:5
    subplot(2,5,i)
    plot( simParsOrdered(:,inds(i,1)),behVarMat(:,inds(i,2)),'.');
    [correls4(1,i,1),correls4(1,i,2)] = corr(simParsOrdered(:,inds(i,1)), ...
        behVarMat(:,inds(i,2)),'type','pearson');
    if p(2) < .05
        l = lsline;
        l.Color = [1 0 0];
        l.LineWidth = 2;
    end
    if i==1
        ylabel(sprintf('%s',varNames{inds(i,2)}),'FontWeight','normal')
    else
        ylabel(varNames{inds(i,2) },'FontWeight','normal')
    end
    if i>2
        xlim([0 1])
    end
    xlabel(parNames{inds(i,1)},'FontWeight','bold')
    title(parNames{inds(i,1)})
    box off
    yl = round(ylim,1, 'significant'); 
    ylim(yl);
    yticks(linspace(yl(1), yl(2), 3));
end

for iPar = 1:4
    subplot(2,5,iPar+6)
    
    plot(simParsOrdered(:,iPar), fitParsOrdered(:,iPar),'.');
    [correls4(2,iPar+1,1),correls4(2,iPar+1,2)] = corr(simParsOrdered(:,iPar), ...
        fitParsOrdered(:,iPar),'type','Spearman');
    if p(2) < .05
        l = lsline;
        l.Color = [1 0 0];
        l.LineWidth = 2;
    end
    if iPar==1
        ylim([0 100])
        ylabel(parNames{iPar},'FontWeight','bold')
        yticks([0 50 100]);
    else
        ylabel(parNames{iPar},'FontWeight','bold')
        if iPar==4
            xlim([0 1])
        end
        ylim([0 1]);
        yticks([0 .5 1]);
    end
    xlabel(parNames{iPar},'FontWeight','bold')
    box off
end

h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(17),'./Figs/MemToolbox2D_behMetrics.jpg')


%% Supplementary figs
clear
load('./Data/MemToolbox1DSim.mat','params2','nll','fitPars','nIter','sd','b','g')

% take some slices across par space

%fix 2 pars, look at correl of other par
legendLab = {'y=x','median','±1SD','±2SD','±3SD','±4SD'};

% guess
figure(18);clf
cols = 3:2:10;
colsSD = cols;
Q = [];
n = length(cols);
parNames = {'\gamma','\beta','\sigma'};
for iB = 1:n
    for iSD=1:n
        subplot(n,n,(iB-1)*n+iSD)
        p = quantPlot(sq(params2(:, 1,cols(iB), cols(iSD),1)), sq(fitPars(:,1,cols(iB),cols(iSD),:)));
        hold on
        line([0 1], [0 1], 'Color','k','LineStyle','--')
        if iB==1
            title(sprintf('\\sigma=%.0f',sd(colsSD(iSD))),'FontWeight','normal')
        end
        if iSD==1
            ylabel(sprintf('\\beta=%.1f',b(cols(iB))))
        end
    end
end
h = findobj(gca, 'Type', 'Line','-or','Type','Patch');
legend(h,legendLab,'Location',[.8064 .1242 .1482 .2012]);

makeSubplotScalesEqual(n,n)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('Simulated \gamma value','FontWeight','normal')
h.YLabel.Visible = 'on';
y = ylabel('Recovered \gamma value','FontWeight','normal');
y.Position(1) = y.Position(1) *1.2;
saveas(figure(18), 'Figs/MemToolbox2D_gamma1D.jpg');

%% misbind
figure(19);clf
for iG = 1:n
    for iSD=1:n
        subplot(n,n,(iG-1)*n+iSD)
        p = quantPlot(sq(params2(cols(iG), 2,:, cols(iSD),1)), sq(fitPars(cols(iG),2,:,cols(iSD),:)));
        hold on
        line([0 1], [0 1], 'Color','k','LineStyle','--')
        if iG==1
            title(sprintf('\\sigma=%.0f',sd(colsSD(iSD))), 'FontWeight','normal')
        end
        if iSD==1
            ylabel(sprintf('\\gamma=%.1f',g(cols(iG))))
        end
        
    end
end
h = findobj(gca, 'Type', 'Line','-or','Type','Patch');
legend(h,legendLab,'Location',[.8064 .1242 .1482 .2012]);

makeSubplotScalesEqual(n,n)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('Simulated \beta value','FontWeight','normal')
h.YLabel.Visible = 'on';
y = ylabel('Recovered \beta value','FontWeight','normal');
y.Position(1) = y.Position(1) *1.2;

saveas(figure(19), 'Figs/MemToolbox2D_beta1D.jpg');
%% SD
figure(20);clf
inds = [];
for iG = 1:n
    for iB=1:n
        if any(~isnan(col(fitPars(cols(iG),3,cols(iB),:,:))))

            subplot(n,n,(iG-1)*n+iB)
            p = quantPlot(sq(params2(cols(iG), 3, cols(iB), :,1)), sq(fitPars(cols(iG),3,cols(iB),:,:)));
            hold on
            line([0 100], [0 100], 'Color','k','LineStyle','--')
            
            if iG==1
                title(sprintf('\\beta=%.1f',b(cols(iB))),'FontWeight','normal')
            end
            if iB==1
                ylabel(sprintf('\\gamma=%.1f',g(cols(iG))))
            end
            axis([0 100 0 100])
            inds = [inds, (iG-1)*n+iB];
        end
    end
end
h = findobj(gca, 'Type', 'Line','-or','Type','Patch');
legend(h,legendLab,'Location',[.7475 .1099 .1482 .2012]);

makeSubplotScalesEqual(n,n,inds)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('Simulated \gamma value','FontWeight','normal')
h.YLabel.Visible = 'on';
y = ylabel('Recovered \gamma value','FontWeight','normal');
y.Position(1) = y.Position(1) *1.2;
saveas(figure(20), 'Figs/MemToolbox2D_sigma1D.jpg');




%%

clear
load('./Data/MemToolbox2DBehMetrics.mat');
varNamesShort = {'TD','NN','TD-NN','SE','SEC','SEm','SECm'};
parNames = {'\sigma','\alpha','\beta','\gamma'};

figure(21)
for iVar = 1:7
    for iPar = 1:4
        subplot(7,4,(iVar-1)*4+iPar)
        plot(simParsOrdered(:,iPar),behVarMat(:,iVar),'.');
        [correls(iVar,iPar,1),correls(iVar,iPar,2)] = corr(simParsOrdered(:,iPar),...
            behVarMat(:,iVar),'type','pearson');
        if correls(iVar,iPar,2) < .05
            l = lsline;
            l.Color = [1 0 0];
            l.LineWidth = 2;
        end
        if iPar==1
            ylabel(varNamesShort{iVar},'FontWeight','normal')
        else
            xlim([0 1])
        end
        if iVar==1
            title(parNames{iPar},'FontWeight','bold')
        elseif iVar==7
            xlabel(parNames{iPar},'FontWeight','bold')
        end
        box off
        y = yticks;
        yticks([y(1) y(end)]);
        x = xticks;
        xticks([x(1) x(end)]);
    end
end

saveas(figure(21), 'Figs/MemToolbox2D_behMetricsAll.jpg');
end