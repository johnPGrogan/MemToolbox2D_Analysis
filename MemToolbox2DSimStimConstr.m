function MemToolbox2DSimStimConstr(nIter)
% MemToolbox2DSimStimConstr(nIter)
% simulation with and within constraints on the minimum distances between stimuli
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
itemsPerTrial   = ones(1,numTrials)*3; % number of items per trial
models          = {SwapModel2D()}; % 1D and 2D misbinding models
nPPs            = 100; % number of pps to simulate
if nargin == 0
    nIter = 100;
end

% distances between [stimuli, stimuli & edges, stimuli & centre]
minDists = [88 29 88]; % pix equiv of [3 3 1] vis deg at 40cm
dimensions = [1366;768]; % screen dimensions 

nSteps = 11;
g = linspace(0.0,1,nSteps);
b = linspace(0,1,nSteps);
sd = linspace(0.01,100,nSteps);

params = [ kron(g, ones(1,nSteps^2)) ;...
    kron(ones(1,nSteps), kron(b, ones(1,nSteps)));...
    kron(ones(1,nSteps^2), sd)]';
params(sum(params(:,1:2),2) >1,:) = [];

params = repmat(params,nIter,1);

nParSets = size(params, 1);
%% fit model

fitPars2D = cell(nPPs,1); %preallocate
ll2D      = fitPars2D; % log likelihoods

tic
parfor i = 1:nParSets
    disp(i)  
    % 2D
    
    % get locations - with min distances
    simulatedData2D{i} = GenerateDisplays2D(numTrials, itemsPerTrial,...
                        1, dimensions); 
    % simulate responses
    simulatedData2D{i}.errors = SampleFromModel2D(models{1}, params(i,:),...
        [1 numTrials], simulatedData2D{i}); % simulate responses
    % fit swap model with maximum likelihood
    [fitPars2D{i,1}, ll2D{i,1}] = MLE(simulatedData2D{i}, models{1});
   
   
    % get locations - with min distances
    simulatedData2D2{i} = GenerateDisplays2D(numTrials, itemsPerTrial,...
                        1, dimensions, minDists);   
    % simulate responses
    simulatedData2D2{i}.errors = SampleFromModel2D(models{1}, params(i,:),...
        [1 numTrials], simulatedData2D2{i}); % simulate responses
    % fit swap model with maximum likelihood
    [fitPars2D2{i,1}, ll2D2{i,1}] = MLE(simulatedData2D2{i}, models{1});
    
end
toc
save('Data/MemToolbox2DSimStimConstr.mat')
%% plot params from swap model

fitPars = cat(3, cat(1,fitPars2D{:}),cat(1, fitPars2D2{:}));

% reorder pars and calcualte a=1-g-B;
fitParsOrdered = [fitPars(:,3,:), 1 - fitPars(:,1,:) - fitPars(:,2,:), fitPars(:,2,:), fitPars(:,1,:)];
simParsOrdered = repmat([params(:,3), 1 - params(:,1) - params(:,2), params(:,2), params(:,1)], 1,1,2);

%%
figure(1);clf
modelNames = {'unconstrained','constrained'};
parNames = {'\sigma','\alpha','\beta','\gamma'}; % standard deviation, target, misbind, guess
for j = 1:2
    for i=1:4
        subplot(2,4,(j-1)*4+i)
		plot(simParsOrdered(:,i,j),fitParsOrdered(:,i,j),'x')
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
%saveas(figure(1), 'Figs/MemToolbox2DSimStimConstr_1.jpg');
%% heatmaps
figure(2);clf

for j = 1:2
    for i=1:4
        subplot(2,4,(j-1)*4+i)
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
figure(3);clf
for j = 1:2
    for i=1:4
        subplot(2,4,(j-1)*4+i)
        s = simParsOrdered(:,i,j);
        
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
%saveas(figure(3), 'Figs/MemToolbox2DSimStimConstr_2.jpg');

%% plot diff in pars vs pars
d = cat(4, fitParsOrdered - simParsOrdered, abs(fitParsOrdered - simParsOrdered));
for i = 1:4
    [~,diffP(i)] = ttest(d(:,i,1), d(:,i,2));
end
labels = {'recovery error','abs( recovery error )'};
figure(4);clf;
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
       title(parNames{i})
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
x = xlabel('simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;

%saveas(figure(4), 'Figs/MemToolbox2DSimStimConstr_3.jpg');
%% bland altman

figure(5);clf
meanVals = (fitParsOrdered + simParsOrdered) ./ 2;
d = fitParsOrdered - simParsOrdered;
c  = get(gca,'ColorOrder');%[1 0 0; 1 .4 0; 1 1 0;0 1 0; 0 1 1; 0 0 1;];

    for i = 1:4
       subplot(1,4,i)
        for j = 1:2
            hold on
            conditionalPlot(meanVals(:,i,j), d(:,i,j),[],'color',c(j,:));
        end
        
       hold on
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
makeSubplotScalesEqual(1,4,[2:4])
legend(h(:,1),modelNames,'Location','Best')
    
%saveas(figure(5), 'Figs/MemToolbox2DSimStimConstr_4.jpg');

%% vs

figure(6);clf
for i = 1:4
    subplot(1,4,i)
    plot(fitParsOrdered(:,i,1), fitParsOrdered(:,i,2), 'x');
    xlabel(modelNames(1))
    if i==1
        ylabel(modelNames(2))
        axis([0 100 0 100])
    else
        axis([0 1 0 1])
    end
    title(parNames{i})
    h = lsline;
    h.Color = [1 0 0];
    h.LineWidth = 2;

end

%%
save('Data/MemToolbox2DSimStimConstr.mat')