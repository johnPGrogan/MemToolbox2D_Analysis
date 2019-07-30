function MemToolbox2DSimNTrials(nIter)
% MemToolbox2DSimNTrials(nIter)
% simulation to see how number of trials affects parameter recovery, in 1D and 2D models
%
% requires the following functions from matlib:
% makeSubplotScalesEqual(), nancat(), conditionalPlot(), sq(),
% 
% nIter is the number of iterations across the parameter sweep (default = 100)
%
% it will save the output in ./Data

%% simulate and look at params
nSteps          = 20;
numTrials       = linspace(10,200,nSteps); % number of trials
% itemsPerTrial   = ones(1,numTrials)*3; % number of items per trial
model           = SwapModel2D(); % 1D and 2D misbinding models
nPPs            = 100; % number of pps to simulate
if nargin == 0
    nIter = 100;
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

fitPars = zeros(nParSets,3,nSteps); %preallocate
fitPars1D = fitPars;
tic
parfor j=1:nSteps
%     fprintf('\n%d ',j)
    disp(j)
    for i = 1:nParSets
        
%         1D
        % get locations
        simulatedData1D = GenerateDisplays(numTrials(j), ones(1,numTrials(j))*3);
        
        % simulate responses
        simulatedData1D.errors = SampleFromModel(SwapModel(), params(i,:),...
            [1 numTrials(j)], simulatedData1D); % simulate responses
        
        % fit swap model with maximum likelihood
        fitPars1D(i,:,j) = MLE(simulatedData1D, SwapModel());
        
        % 2D
        
        % get locations
        simulatedData2D = GenerateDisplays2D(numTrials(j), ones(1,numTrials(j))*3);
        
        % simulate responses
        simulatedData2D.errors = SampleFromModel2D(model, params(i,:),...
            [1 numTrials(j)], simulatedData2D); % simulate responses
        
        % fit swap model with maximum likelihood
        fitPars(i,:,j)= MLE(simulatedData2D, model);
        
        pars{i,j} = params(i,:);
        
    end
end
toc
save('Data/MemToolbox2DSimNTrials.mat')

%% plot params from swap model

for i = 1:3
    simPars(:,i,:) = cellfun(@(x) x(:,i), pars);
end

% reorder pars and calcualte a=1-g-B;
fitParsOrdered = [fitPars(:,3,:), 1 - fitPars(:,1,:) - fitPars(:,2,:), fitPars(:,2,:), fitPars(:,1,:)];
simParsOrdered = repmat([simPars(:,3), 1 - simPars(:,1) - simPars(:,2), simPars(:,2), simPars(:,1)], 1,1,1);

fitPars1DOrdered = [fitPars1D(:,3,:), 1 - fitPars1D(:,1,:) - fitPars1D(:,2,:), fitPars1D(:,2,:), fitPars1D(:,1,:)];

%%
figure()

m = ['ko';'bx';'g+';'r^'];
modelNames = {'1D','2D'};
parNames = {'SD','\alpha','\beta','\gamma'}; % standard deviation, target, misbind, guess
for j = 1:nSteps
    subplot(4,5,j)
    hold on
    plot(simParsOrdered(:,1,1)./100,fitParsOrdered(:,1,j)./100,m(1,:))
    for k=1:3
        hold on
        plot(simParsOrdered(:,k+1,1),fitParsOrdered(:,k+1,j),m(k+1,:))
    end
    lsline
    title(num2str(numTrials(j)))
    
end

%% abs mean error

absMeanErr = abs(nancat(4, fitParsOrdered, fitPars1DOrdered) - simParsOrdered) ./ simParsOrdered;
titles = {'1D', '2D'};
figure()
c = get(gca,'ColorOrder');

for k = 1:2
    for j = 1:2
        subplot(2,2,(j-1)*2+k) 
        if j==1
            a = absMeanErr(:,:,:,k);
        else % take good pars
            a = absMeanErr(137:726:nParSets,:,:,k);
        end

        h = errorBarPlot(permute(a(:,2:4,:),[1,3,2]),'area',1);
        hold on
        yyaxis right
        h2 = errorBarPlot(permute(a(:,1,:),[1,3,2]),'area',1);
        gc = gca;
        gc.YColor = c(1,:);
        h2(1,1).Color = c(1,:);
        h2(1,1).LineStyle = '-';
        h2(1,2).FaceColor = c(1,:);
        for i = 1:3
            h(i,1).Color = c(i+1,:);
            h(i,2).FaceColor = c(i+1,:);
        end
        if i==2, ylabel('SD abs mean err'), end
        yyaxis left
        gc = gca;
        gc.YColor = [0 0 0];
        
        xlim([1 20])
        set(gca,'XTick',1:4:20,'XTickLabel',10:(nSteps*2):200)
        if i==1, ylabel('abs mean error'), end
        if j==2, xlabel('number of trials'), end
        if j==1, title(titles{k}), end
    end
end
legend([h2(:,1);h(:,1)],parNames,'Location','NorthEast')
%saveas(figure(2), 'Figs/MemToolbox2DSimNTrials_1.jpg');
%%
figure(3);clf
cols = [1 2 4 8 16 20];
for j = 1:length(cols)
    for i=1:4
        subplot(length(cols),4,(j-1)*4+i)
        s = simParsOrdered(:,i);
        
        f = fitParsOrdered(:,i,cols(j));
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
            ylabel(numTrials(cols(j)))
            if j==1
                axis([0 30 0 30])
            else
                axis([0 100 0 100])
            end
        else
            axis([0 1 0 1])
        end
    end
end
%saveas(figure(3), 'Figs/MemToolbox2DSimNTrials_2.jpg');

%% plot diff in pars vs pars
d = nancat(4, fitPars1DOrdered - simParsOrdered, fitParsOrdered - simParsOrdered);
parNames = {'\sigma','\alpha','\beta','\gamma'};
% a = 1 - g - b;
figure(4);clf;
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
xlabel('simulated parameter','FontWeight','bold')

%saveas(figure(4), 'Figs/MemToolbox2DSimNTrials_3.jpg');

%% bland altman

figure(5);clf
meanVals = nancat(4, fitPars1DOrdered + simParsOrdered, fitParsOrdered + simParsOrdered) ./ 2;

for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
       set(gca,'ColorOrder',c);
        for k = 1:20
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
           xlim([0 1]);
       end
       box off
    end
end
makeSubplotScalesEqual(2,4,[2:4, 6:8])
makeSubplotScalesEqual(2,4,[1 5])
colormap(c)
subplot(2,4,8), h1 = colorbar(gca,'North','Ticks',[0 1], 'TickLabels',[numTrials([1 end])]);
    
%saveas(figure(5), 'Figs/MemToolbox2DSimNTrials_4.jpg');

%% plot diff in pars vs pars
figure(6);clf;
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
xlabel('simulated parameter','FontWeight','normal')

%saveas(figure(6), 'Figs/MemToolbox2DSimNTrials_5.jpg');
%%
save('./Data/MemToolbox2DSimNTrials.mat')