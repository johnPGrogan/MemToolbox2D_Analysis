function MemToolbox2DSim2AFC(nIter)
% MemToolbox2DSim2AFC(nIter)
% simulation to examine 2AFC version of model. Uses the StandardMixtureModel for speed
%
% requires the following functions from matlib:
% makeSubplotScalesEqual(), nancat(), conditionalPlot(), sq(),
% 
% nIter is the number of iterations across the parameter sweep (default = 100)
%
% it will save the output in ./Data
%


%% simulate and look at params

numTrials       = 100; % number of trials
itemsPerTrial   = ones(1,numTrials); % number of items per trial
models          = {TwoAFC(StandardMixtureModel()), TwoAFC2D(StandardMixtureModel2D())}; % models
nPPs            = 100; % number of pps to simulate
if nargin == 0
    nIter = 100;
end

nSteps = 11;
g = linspace(0.0,1,nSteps);
sd = linspace(0.01,100,nSteps);

params = [ kron(g, ones(1,nSteps));...
    kron(ones(1,nSteps), sd)]';

params = repmat(params,nIter,1);

nParSets = size(params, 1);
%% fit model

fitPars1D = cell(nPPs,1); %preallocate
ll1D      = cell(nPPs,1); % log likelihoods
fitPars2D = fitPars1D; %preallocate
ll2D      = ll1D; % log likelihoods

tic
parfor i = 1:nParSets
    disp(i)
    % 1D
    
    % get locations
    simulatedData1D{i} = GenerateDisplays(numTrials, itemsPerTrial, 2); 
    
    % simulate responses
    simulatedData1D{i}.afcCorrect = SampleFromModel(models{1}, params(i,:),...
        [1 numTrials], simulatedData1D{i})'; % simulate responses
    
    % fit swap model with maximum likelihood
    [fitPars1D{i,1}, ll1D{i,1}] = MLE(simulatedData1D{i}, models{1});
    
    
    % 2D
    
    % get locations
    simulatedData2D{i} = GenerateDisplays2D(numTrials, itemsPerTrial, 2, [1366;768]); 
    
    % simulate responses
    simulatedData2D{i}.afcCorrect = SampleFromModel2D(models{2}, params(i,:),...
        [1 numTrials], simulatedData2D{i}); % simulate responses
    
    % fit swap model with maximum likelihood
    [fitPars2D{i,1}, ll2D{i,1}] = MLE(simulatedData2D{i}, models{2});
   
    
end
t = toc
save('Data/MemToolbox2DSim2AFC.mat')
%% plot params from swap model

fitPars = cat(3, cat(1,fitPars1D{:}), cat(1,fitPars2D{:}));

% reorder pars and calcualte a=1-g-B;
fitParsOrdered = [fitPars(:,2,:), 1 - fitPars(:,1,:),  fitPars(:,1,:)];
simParsOrdered = repmat([params(:,2), 1 - params(:,1), params(:,1)], 1,1,2);

%%
modelNames = {'1D','2D'};
parNames = {'SD','\alpha','\gamma'};
for j = 1:2
    for i=1:3
        subplot(2,3,(j-1)*3+i)
        if j==2 && i==1
            plot(simParsOrdered(:,i,j),fitParsOrdered(:,i,j),'x')
        else
            plot(simParsOrdered(:,i,j),fitParsOrdered(:,i,j),'x')
        end
        lsline
        if j==1
            title(parNames{i})
        else
            xlabel('sim pars')
        end
        if i==1
            ylabel(modelNames{j})
            axis([0 100 0 100])
        else
            axis([0 1 0 1])
        end
    end
end
%saveas(figure(1), 'Figs/MemToolbox2DSim2AFC_1.jpg');
%% heatmaps
figure()

for j = 1:2
    for i=1:3
        subplot(2,3,(j-1)*3+i)
        h = heatmap(simParsOrdered(:,i,j), fitParsOrdered(:,i,j));
        if j==1
            title(parNames{i})
        else
            xlabel('sim pars')
        end
        h.YData = h.YData([2,1]);
        gc = gca;
        set(gca,'YTick',gc.XTick-1,'YTickLabel',gc.XTickLabel(end:-1:1))
        if i==1
            ylabel(modelNames{j})
            gc.XTickLabel = cellfun(@(x) {str2num(x)},gc.XTickLabel);
            set(gca,'XTickLabel',gc.XTickLabel)
        else
%             axis([0 1 0 1])
        end

    end
end

%%
figure(3)
for j = 1:2
    for i=1:3
        subplot(2,3,(j-1)*3+i)
        if j==2 && i==1
             s = simParsOrdered(:,i,j);
        else
            s = simParsOrdered(:,i,j);
        end
        
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
            ylabel(modelNames{j})
			axis([0 100 0 100])
        else
            axis([0 1 0 1])
        end
    end
end
%saveas(figure(3), 'Figs/MemToolbox2DSim2AFC_2.jpg');


%% plot diff in pars vs pars
d = fitParsOrdered - simParsOrdered;

figure(4);clf;
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
x = xlabel('simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2) .* 1.2;
%saveas(figure(4), 'Figs/MemToolbox2DSim2AFC_3.jpg');

%% bland altman

figure(5);clf
meanVals = (fitParsOrdered + simParsOrdered) ./ 2;

for j = 1:2
    for i = 1:3
       subplot(2,3,(j-1)*3+i)
        for k = 1
            conditionalPlot(meanVals(:,i,j), d(:,i,j),[]);%,'color',c(k,:));
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
makeSubplotScalesEqual(2,3,[2:3, 5:6])
makeSubplotScalesEqual(2,3,[1 4])
    
%saveas(figure(5), 'Figs/MemToolbox2DSim2AFC_4.jpg');

%%
save('Data/MemToolbox2DSim2AFC.mat')