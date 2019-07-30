function MemToolbox2DSimConstBias(nIter)
% MemToolbox2DSimConstBias(nIter)
% simulation to see how a constant translation bias affects the model fits, 
% and how the bias correction functions works for it
%
% requires the following functions from matlib:
% makeSubplotScalesEqual(), nancat(), conditionalPlot(), sq(),
% 
% nIter is the number of iterations across the parameter sweep (default = 100)
%
% it will save the output in ./Data
%% simulate model

numTrials       = 100; % number of trials
itemsPerTrial   = ones(1,numTrials)*3; % number of items per trial
model           = SwapModel2D(); % 2D misbinding model
nPPs            = 100; % number of pps to simulate
if nargin == 0
    nIter = 100;
end
nSteps = 6; % number of steps in bias size
bias = linspace(20,200,nSteps); % different bias sizes to simulate

%% grid search of params
nParSteps = 11;
cols = [3:2:10];
g = linspace(0.1,.98,nParSteps);
b = linspace(0.1,.98,nParSteps);
sd = linspace(0.01,100,nParSteps);

params = [ kron(g, ones(1,nParSteps^2)) ;...
    kron(ones(1,nParSteps), kron(b, ones(1,nParSteps)));...
    kron(ones(1,nParSteps^2), sd)]';

params(sum(params(:,1:2),2) >1,:) = [];

params = repmat(params,nIter,1);

nParSets = length(params);
dimensions = [1366; 768]; % screen dims
%%
[fitPars1, fitPars2,fitPars3] = deal(NaN(nParSets,3,nSteps));
fitPars3(:,4:5,:) = NaN;

tic
parfor j=1:nSteps
    disp(j)

        for i=1:nParSets
            % get locations
            simData = GenerateDisplays2D(numTrials, itemsPerTrial);

            % simulate responses
            simData.errors = SampleFromModel2D(model, params(i,:),...
                [1 numTrials], simData); % simulate responses% simulate responses

            % fit swap model with maximum likelihood
            fitPars1(i,:,j) = MLE(simData, model);


            % shift resps towards bottom left hand side
            resps = simData.errors + simData.targets;
            resps2 = resps - bias(j); % make shift towards bototm left
            if any(resps2 < 0,'all')
                resps2(resps2<0) = 0;
            end
            err2 = resps2 - simData.targets; % new errors
            
            simData2 = simData;
            simData2.errors = err2;
            fitPars2(i,:,j) = MLE(simData2, model);
            
            % correct for bias
            fitPars3(i,:,j) = MLE(simData2, WithBias2D(model));
        end
        
%     end
end
t = toc
save('Data/MemToolbox2DSimConstBias.mat')
%%

simPars = params;

% reorder pars and calcualte a=1-g-B;
fitParsOrdered = [fitPars2(:,3,:), 1 - fitPars2(:,1,:) - fitPars2(:,2,:), fitPars2(:,2,:), fitPars2(:,1,:)];
simParsOrdered = repmat([simPars(:,3), 1 - simPars(:,1) - simPars(:,2), simPars(:,2), simPars(:,1)], 1,1,nSteps);


fitParsOrdered3 =  [fitPars3(:,3,:), 1 - fitPars3(:,1,:) - fitPars3(:,2,:), fitPars3(:,2,:), fitPars3(:,1,:), fitPars3(:,4,:),fitPars3(:,5,:)];

%%
figure()

m = ['ko';'bx';'g+';'r^'];
modelNames = {'1D','2D'};
parNames = {'SD','\alpha','\beta','\gamma'};
for j = 1:nSteps
    subplot(3,2,j)
    hold on
    plot(simParsOrdered(:,1,1)./100,fitParsOrdered(:,1,j)./100,m(1,:))
    for k=1:3
        hold on
        plot(simParsOrdered(:,k+1,1),fitParsOrdered(:,k+1,j),m(k+1,:))
    end
    lsline
    title(num2str(bias(j)))
    
end

%% abs mean error

absMeanErr{1} = abs(fitParsOrdered - simParsOrdered) ./ simParsOrdered;
absMeanErr{2} = abs(fitParsOrdered3(:,1:4,:) - simParsOrdered) ./ simParsOrdered;
labels = {'bias','corrected'};
figure(2);clf
c = get(gca,'ColorOrder');
for j=1:2
    subplot(1,2,j)
    h = errorBarPlot(permute(absMeanErr{j}(:,2:4,:),[1,3,2]),'area',1);
    hold on
    yyaxis right
    h2 = errorBarPlot(permute(absMeanErr{j}(:,1,:),[1,3,2]),'area',1);
    gc = gca;
    gc.YColor = c(1,:);
    h2(1,1).Color = c(1,:);
    h2(1,1).LineStyle = '-';
    h2(1,2).FaceColor = c(1,:);
    for i = 1:3
        h(i,1).Color = c(i+1,:);
        h(i,2).FaceColor = c(i+1,:);
    end
    ylim([0 6000])
    yyaxis left
    gc = gca;
    gc.YColor = [0 0 0];

    xlim([1 nSteps])
    set(gca,'XTick',1:nSteps,'XTickLabel',bias)
    ylabel('abs mean error')
    xlabel('bias')
    
    title(labels{j})
    
end
legend([h2(:,1);h(:,1)],parNames,'Location','Best')

makeSubplotScalesEqual(1,2)
%saveas(figure(2), 'Figs/MemToolbox2DSimConstBias_1.jpg');
%%
figure(3);clf
for j = 1:nSteps
    for i=1:4
        subplot(nSteps,4,(j-1)*4+i)
        s = simParsOrdered(:,i);        
        f = fitParsOrdered(:,i,j);
        fp = groupMeans(f,1,s,'dim');
        quantPlot(unique(s),fp);
        hold on
        line([0 100],[0 100],'Color','k','LineStyle','--')
        if j==1
            title(parNames{i})
        else
            xlabel('sim pars')
        end
        if i==1
            ylabel(num2str(bias(j)))
            axis([0 100 0 100])
        else
            axis([0 1 0 1])
        end
    end
end
SuperTitle('bias')
%saveas(figure(3), 'Figs/MemToolbox2DSimConstBias_2.jpg');

%%
figure(4);clf
for j = 1:nSteps
    for i=1:4
        subplot(nSteps,4,(j-1)*4+i)
        s = simParsOrdered(:,i);        
        f = fitParsOrdered3(:,i,j);
        fp = groupMeans(f,1,s,'dim');
        quantPlot(unique(s),fp);
        hold on
        line([0 100],[0 100],'Color','k','LineStyle','--')
        if j==1
            title(parNames{i})
        else
            xlabel('sim pars')
        end
        if i==1
            ylabel(num2str(bias(j)))
            axis([0 100 0 100])
        else
            axis([0 1 0 1])
        end
    end
end
SuperTitle('corrected')
%saveas(figure(4), 'Figs/MemToolbox2DSimConstBias_3.jpg');

%% plot diff in pars vs pars
d = nancat(4, fitParsOrdered - simParsOrdered, fitParsOrdered3(:,1:4,:) - simParsOrdered);

figure(5);clf;
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
xlabel('simulated parameter','FontWeight','normal')
%saveas(figure(5), 'Figs/MemToolbox2DSimConstBias_4.jpg');

%% bland altman

figure(6);clf
meanVals = nancat(4, fitParsOrdered + simParsOrdered, fitParsOrdered3(:,1:4,:) + simParsOrdered) ./ 2;

for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
        for k = 1:nSteps
            conditionalPlot(meanVals(:,i,k,j), d(:,i,k,j),[],'color',c(k,:));
            hold on; 
        end
        line([0 100],[0 0],'Color','k','LineStyle','--')

       if j==1; title(parNames{i}); end
       if j==2,    xlabel('mean (fit + sim)'), end
       if i==1
           ylabel(sprintf('%s\nfit - sim',modelNames{j}));
           xlim([0 100]);
       else
           xlim([0 1])
       end
       
       box off
    end
end
makeSubplotScalesEqual(2,4,[2:4, 6:8])
makeSubplotScalesEqual(2,4,[1 5])
c = c(1:nSteps,:);
colormap(c)
subplot(2,4,8), colorbar(gca,'South','Ticks',[0 1], 'TickLabels',arrayfun(@num2str, bias([1 end]),'UniformOutput',0))
    
%saveas(figure(6), 'Figs/MemToolbox2DSimConstBias_5.jpg');

%% plot diff in pars vs pars
d = abs(nancat(4, fitParsOrdered - simParsOrdered, fitParsOrdered3(:,1:4,:) - simParsOrdered));

figure(7);clf;
simParsOrdered1 = round(simParsOrdered,2,'significant');
c  = [ 0 0 1; 0 1 1;0 1 0; 1 1 0; 1 .4 0;1 0 0;];
c = c(1:nSteps,:);
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
xlabel('simulated parameter','FontWeight','normal')
%saveas(figure(7), 'Figs/MemToolbox2DSimConstBias_6.jpg');
%%

save('Data/MemToolbox2DSimConstBias.mat')