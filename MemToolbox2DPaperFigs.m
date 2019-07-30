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

saveas(figure(1), 'Figs/MemToolbox2D_1.jpg');
%% misbind
figure(2);clf
inds = [];
for iG = 1:n
    for iSD=1:n
        if any(~isnan(col(fitPars(cols(iG),2,:,cols(iSD),:))))
            subplot(n,n,(iG-1)*n+iSD)
%             h = heatmap(col(params2(cols(iG),2,:,cols(iSD),:)), col(fitPars(cols(iG),2,:,cols(iSD),:)));
%             plot(col(params2(cols(iG),2,:,cols(iSD),:)), col(fitPars(cols(iG),2,:,cols(iSD),:)),'x')
            p = quantPlot(sq(params2(cols(iG), 2,:, cols(iSD),1)), sq(fitPars(cols(iG),2,:,cols(iSD),:)));
            hold on
            line([0 1],[0 1],'Color','k','LineStyle','--')
            if iG==1
                title(sprintf('\\sigma=%.1f',sd(cols(iSD))),'FontWeight','normal')
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
xlabel('Simulated \beta value','FontWeight','normal')
h.YLabel.Visible = 'on';
y = ylabel('Recovered \beta value','FontWeight','normal');
y.Position(1) = y.Position(1) *1.2;

saveas(figure(2), 'Figs/MemToolbox2D_2.jpg');
%% SD
figure(3);clf
inds = [];
for iG = 1:n
    for iB=1:n
        if any(~isnan(col(fitPars(cols(iG),3,cols(iB),:,:))))
            subplot(n,n,(iG-1)*n+iB)
%             h = heatmap(col(params2(cols(iG),3,cols(iB),:,:)), col(fitPars(cols(iG),3,cols(iB),:,:)));
%             plot(col(params2(cols(iG),3,cols(iB),:,:)), col(fitPars(cols(iG),3,cols(iB),:,:)),'x')
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
xlabel('Simulated \sigma value','FontWeight','normal')
h.YLabel.Visible = 'on';
y = ylabel('Recovered \sigma value','FontWeight','normal');
y.Position(1) = y.Position(1) *1.2;

saveas(figure(3), 'Figs/MemToolbox2D_3.jpg');


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
       x = unique(round(simParsOrdered1(:,i,j), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
    end
end
makeSubplotScalesEqual(2,4,[2:4, 6:8])
makeSubplotScalesEqual(2,4,[1 5])
colormap(c)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(4), 'Figs/MemToolbox2D_4.jpg');

%%
clear
load('./Data/MemToolbox2DPosteriors.mat', 'mcmc')
modelNames = {'1D';'2D'};
set(0,'DefaultAxesFontWeight','bold')
f = figure(5);
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
        if j==1
            title([parNames{inds(i,1)} ' & ' parNames{inds(i,2)}])
        end

        [correls(j,i,1),correls(j,i,2)] = corr(mcmc(j).vals(:,inds(i,1)), mcmc(j).vals(:,inds(i,2)), 'type', 'Spearman');
    end
end
colormap(palettablecolormap('sequential'));
makepalettable(f);

disp(correls)
saveas(figure(5),'Figs/MemToolbox2D_5.jpg')

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
       x = unique(round(simParsOrdered1(:,i), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
    end
end
subplot(2,4,8), colorbar(gca,'South','Ticks',[0 1], 'TickLabels',[numTrials([1 end])])

makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])

h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(6), 'Figs/MemToolbox2D_6.jpg');

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
       if i==1, ylabel(sprintf('%s\nAbs( recovery error )',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
    end
end
subplot(2,4,8), colorbar(gca,'North','Ticks',[0 1], 'TickLabels',[numTrials([1 end])])


makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(7), 'Figs/MemToolbox2D_7.jpg');


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
       if j==1; title(parNames{i}); end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(sprintf('%s\nrecovery error',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
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

saveas(figure(8), 'Figs/MemToolbox2D_8.jpg');

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
       if i==1, ylabel(sprintf('%s\nrecovery error',labels{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x([1 ceil(length(x)/2) end]))
       xlim([1 11])
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

saveas(figure(9), 'Figs/MemToolbox2D_9.jpg');

%%

clear
load('./Data/MemToolbox2DSimRadBias.mat','fitParsOrdered','fitPars3Ordered',...
    'simParsOrdered','nSteps','bias','nParSets','modelNames','parNames')
labels = {'Uncorrected','Corrected'};

d = nancat(5, fitParsOrdered(:,1:4,:) - simParsOrdered, fitPars3Ordered(:,1:4,:) - simParsOrdered);
parNames = {'\sigma','\alpha','\beta','\gamma'};


figure(10);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c = colormap('jet');
c = c(1:11:end,:);
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
           ylabel(sprintf('%s\nrecovery error',labels{j}),'FontWeight','bold')
       end
    end
end
subplot(2,4,5), colorbar(gca,'North','Ticks',[0 .5 1],'TickLabels',[bias(1) 0 bias(end)])
makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(10), 'Figs/MemToolbox2D_10.jpg');


%%

clear

load('./Data/MemToolbox2DSimDistrRadial.mat','fitParsOrdered','fitPars3Ordered','simParsOrdered','bias','nIter')
labels = {'Uncorrected','Corrected'};
% plot diff in pars vs pars
d = nancat(5, fitParsOrdered(:,1:3,:) - simParsOrdered, fitPars3Ordered(:,1:3,:) - simParsOrdered);
parNames = {'\sigma','\alpha','\gamma'};


figure(11);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c = colormap('jet');
c = c(1:11:end,:);
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
       box off
       x = unique(round(simParsOrdered1(:,i), 1));
       set(gca,'XTick',1:4:11,'XTickLabel',x([1 ceil(length(x)/2) end]))
       xlim([1 11])
       if i==1
           ylabel(sprintf('%s\nrecovery error',labels{j}),'FontWeight','bold')
       end
    end
end
subplot(2,3,4), colorbar(gca,'North','Ticks',[0 .5 1],'TickLabels',[bias(1) 0 bias(end)])
makeSubplotScalesEqual(2,3,[2:3 5:6])
makeSubplotScalesEqual(2,3,[1 4])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2) .* 1.2;

saveas(figure(11), 'Figs/MemToolbox2D_11.jpg');

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
figure(12);clf;
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
       if j==1, title(parNames{i}), end
       xlabel(parNames{i},'FontWeight','bold')
       if i==1, ylabel(labels{j},'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i,1), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
    end
end
makeSubplotScalesEqual(2,4,[2:4])
makeSubplotScalesEqual(2,4,[6:8])
legend(h(:,1),modelNames,'Location','Best')
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;
saveas(figure(12), 'Figs/MemToolbox2D_12.jpg');

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
       if j==1, title(parNames{i}); end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(labels{j},'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i,1), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
    end
end
makeSubplotScalesEqual(2,4,[2:4])
makeSubplotScalesEqual(2,4,[6:8])
legend(h(:,1),modelNames,'Location','Best')
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;
saveas(figure(13), 'Figs/MemToolbox2D_13.jpg');

%%
clear
load('./Data/MemToolbox2DSim2AFC.mat','fitParsOrdered','simParsOrdered','modelNames')
parNames = {'\sigma','\alpha','\gamma'};
% plot diff in pars vs pars
d = fitParsOrdered - simParsOrdered;

figure(14);clf;
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
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(sprintf('%s\nrecovery error',modelNames{j}),'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i), 1));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
    end
end
makeSubplotScalesEqual(2,3,[2:3, 5:6])
makeSubplotScalesEqual(2,3,[1 4])
colormap(c)
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2) .* 1.2;

saveas(figure(14), 'Figs/MemToolbox2D_14.jpg');


%%
clear
load('./Data/MemToolbox2DFallonPosteriors.mat');

modelNames = {'1D';'2D'};
f = figure(15);clf;
inds = [1 2; 1 3; 2 3;];
parNames = {'\gamma', '\beta','\sigma'};
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
            p(k,i,j) = plot(NaN, NaN,'s','Color',colours(k),...
                'MarkerFaceColor',colours(k),'MarkerSize',4); % plot NaN so colours show up in legend
        end
        
        box off
        axis tight;
        box off
        if j==1
            title([parNames{inds(i,1)} ' & ' parNames{inds(i,2)}],'FontWeight','bold');
        else
            xlabel(parNames{inds(i,1)},'FontWeight','bold')
        end
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
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2) .* 1.2;

saveas(figure(15),'./Figs/MemToolbox2D_15.jpg')


%%

clear
load('./Data/MemToolbox2DBehMetrics.mat');
varNames = {'targDist', 'nearestNeighbour','targDist - nearestNeighbour','swapError','swapErrorCorr','swapErrorMeanThresh','swapErrorMeanThreshCorr'};
parNames = {'\sigma','\alpha','\beta','\gamma'};

figure(16);
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
    else
        ylabel(parNames{iPar},'FontWeight','bold')
        if iPar==4
            xlim([0 1])
        end
    end
    xlabel(parNames{iPar},'FontWeight','bold')
    

    box off
end

h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

saveas(figure(16),'./Figs/MemToolbox2D_16.jpg')


%% Supplementary figs
clear
load('./Data/MemToolbox1DSim.mat','params2','nll','fitPars','nIter','sd','b','g')

% take some slices across par space

%fix 2 pars, look at correl of other par
legendLab = {'y=x','median','±1SD','±2SD','±3SD','±4SD'};

% guess
figure(17);clf
cols = 3:2:10;
colsSD = cols;
Q = [];
n = length(cols);
parNames = {'g','b','sd'};
for iB = 1:n
    for iSD=1:n
        subplot(n,n,(iB-1)*n+iSD)
        p = quantPlot(sq(params2(:, 1,cols(iB), cols(iSD),1)), sq(fitPars(:,1,cols(iB),cols(iSD),:)));
        hold on
        line([0 1], [0 1], 'Color','k','LineStyle','--')
        if iB==1
            title(sprintf('\\sigma=%.1f',sd(colsSD(iSD))),'FontWeight','normal')
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
saveas(figure(17), 'Figs/MemToolbox2D_17.jpg');
%% misbind
figure(18);clf
for iG = 1:n
    for iSD=1:n
        subplot(n,n,(iG-1)*n+iSD)
        p = quantPlot(sq(params2(cols(iG), 2,:, cols(iSD),1)), sq(fitPars(cols(iG),2,:,cols(iSD),:)));
        hold on
        line([0 1], [0 1], 'Color','k','LineStyle','--')
        if iG==1
            title(sprintf('\\sigma=%.1f',sd(colsSD(iSD))), 'FontWeight','normal')
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

saveas(figure(18), 'Figs/MemToolbox2D_18.jpg');
%% SD
figure(19);clf
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
saveas(figure(19), 'Figs/MemToolbox2D_19.jpg');




%%

clear
load('./Data/MemToolbox2DBehMetrics.mat');
varNamesShort = {'TD','NN','TD-NN','SE','SEC','SEm','SECm'};
parNames = {'\sigma','\alpha','\beta','\gamma'};

figure(20)
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
    end
end

saveas(figure(20), 'Figs/MemToolbox2D_20.jpg');
end