function MemToolbox2DSimNDist(nIter)
% MemToolbox2DSimNDist(nIter)
% simulation to see how the number of distractors affects parameter recovery in 1D and 2D models
%
% requires the following functions from matlib:
% makeSubplotScalesEqual(), nancat(), conditionalPlot(), sq(),
% 
% nIter is the number of iterations across the parameter sweep (default = 100)
%
% it will save the output in ./Data
%

%% simulate and look at params

nMaxDist        = 6; % max number of distractors
numTrials       = 100; % number of trials
model           = SwapModel2D(); % 1D and 2D misbinding models
if nargin==0
    nIter        = 100; % number of pps to simulate
end


%% generate random pars with constraint that g+B <=1
nParSteps = 11;
g = linspace(0.01,0.98,nParSteps);
b = linspace(0.01,0.98,nParSteps);
sd = linspace(0.01,100,nParSteps);

params = [ kron(g, ones(1,nParSteps^2)) ;...
    kron(ones(1,nParSteps), kron(b, ones(1,nParSteps)));...
    kron(ones(1,nParSteps^2), sd)]';
params(sum(params(:,1:2),2) >1,:) = [];

params = repmat(params,nIter,1);

nParSets = size(params, 1);


%% fit model

fitPars = zeros(nParSets,3,nMaxDist); %preallocate
fitPars1D = fitPars;
tic
parfor j=1:nMaxDist
    disp(j)
    for i = 1:nParSets
%           1D
        % get locations
        simulatedData1D = GenerateDisplays(numTrials, ones(1,numTrials)*(j+1));
        
        % simulate responses
        simulatedData1D.errors = SampleFromModel(SwapModel(), params(i,:),...
            [1 numTrials], simulatedData1D); % simulate responses
        
        % fit swap model with maximum likelihood
        fitPars1D(i,:,j) = MLE(simulatedData1D, SwapModel());
        
        % 2D   
        % get locations
        simulatedData2D = GenerateDisplays2D(numTrials, ones(1,numTrials)*(j+1));
        
        % simulate responses
        simulatedData2D.errors = SampleFromModel2D(model, params(i,:),...
            [1 numTrials], simulatedData2D); % simulate responses
        
        % fit swap model with maximum likelihood
        [fitPars(i,:,j), ll(i,j)] = MLE(simulatedData2D, model);
        
        pars{i,j} = params(i,:);
    end
end
t=toc
save('Data/MemToolbox2DSimNDist.mat');
%% plot params from swap model
for i = 1:3
    simPars(:,i,:) = cellfun(@(x) x(:,i), pars);
end

% reorder pars and calcualte a=1-g-B;
fitParsOrdered = [fitPars(:,3,:), 1 - fitPars(:,1,:) - fitPars(:,2,:), fitPars(:,2,:), fitPars(:,1,:)];
simParsOrdered = repmat([simPars(:,3), 1 - simPars(:,1) - simPars(:,2), simPars(:,2), simPars(:,1)], 1,1,1);

fitPars1DOrdered = [fitPars1D(:,3,:), 1 - fitPars1D(:,1,:) - fitPars1D(:,2,:), fitPars1D(:,2,:), fitPars1D(:,1,:)];


%%
figure(1);clf

m = ['ko';'bx';'g+';'r^'];
modelNames = {'1D','2D'};
parNames = {'SD','\alpha','\beta','\gamma'}; % standard deviation, target, misbind, guess
for j = 1:nMaxDist
    for i =1:4
        subplot(nMaxDist,4,(j-1)*4+i)
        if i==1
            s = round(simParsOrdered(:,i)./10);
        else
            s = round(simParsOrdered(:,i),1);
        end
        
        f = fitParsOrdered(:,i,j);
        fp = groupMeans(f,1,s,'dim');
        if i==1
            s=s.*10;
        end
        quantPlot(unique(s),fp);
        hold on
        line([0 100],[0 100],'Color','k','LineStyle','--')
        if j==1
            title(parNames{i})
        end
        if i==1
            ylabel(num2str(j))
            axis([0 100 0 100])
        else
            axis([0 1 0 1])
        end
    end
end


%saveas(figure(1), 'Figs/MemToolbox2DSimNDist_1.jpg');


%% 1D

figure(2);clf

m = ['ko';'bx';'g+';'r^'];
modelNames = {'1D','2D'};
parNames = {'SD','\alpha','\beta','\gamma'}; % standard deviation, target, misbind, guess
for j = 1:nMaxDist
    for i =1:4
        subplot(nMaxDist,4,(j-1)*4+i)
        if i==1
            s = round(simParsOrdered(:,i)./10);
        else
            s = round(simParsOrdered(:,i),1);
        end
        
        f = fitPars1DOrdered(:,i,j);
        fp = groupMeans(f,1,s,'dim');
        if i==1
            s=s.*10;
        end
        quantPlot(unique(s),fp);
        hold on
        line([0 100],[0 100],'Color','k','LineStyle','--')
        if j==1
            title(parNames{i})
        end
        if i==1
            ylabel(num2str(j))
            axis([0 100 0 100])
        else
            axis([0 1 0 1])
        end
    end
end
SuperTitle('1D')

%saveas(figure(2), 'Figs/MemToolbox2DSimNDist_2.jpg');

%% plot diff in pars vs pars
modelNames = {'1D', '2D'};
d = nancat(4, fitPars1DOrdered - simParsOrdered, fitParsOrdered - simParsOrdered);

figure(4);clf;
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
       title(parNames{i})
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
xlabel('simulated parameter')


%saveas(figure(4), 'Figs/MemToolbox2DSimNDist_4.jpg');

%% bland altman

figure(5);clf
meanVals = nancat(4, fitPars1DOrdered + simParsOrdered, fitParsOrdered + simParsOrdered) ./ 2;
simParsOrdered1 = round(meanVals,2,'significant');

for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
        for k = 1:6
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
colormap(c)
subplot(2,4,8), h1 = colorbar(gca,'North','Ticks',[0 1], 'TickLabels',[1 nMaxDist]);
    
%saveas(figure(5), 'Figs/MemToolbox2DSimNDist_5.jpg');

%% plot diff in pars vs pars
d = abs(nancat(4, fitPars1DOrdered - simParsOrdered, fitParsOrdered - simParsOrdered));

figure(6);clf;
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
       title(parNames{i})
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
subplot(2,4,8), colorbar(gca,'North','Ticks',[0 1], 'TickLabels',[1 nMaxDist])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('simulated parameter')

%saveas(figure(6), 'Figs/MemToolbox2DSimNDist_6.jpg');
%%
save('Data/MemToolbox2DSimNDist.mat')