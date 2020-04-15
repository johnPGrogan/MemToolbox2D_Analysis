function MemToolbox2DSimPropBias(nIter)
% MemToolbox2DSimPropBias(nIter)
% simulation with radial response biases towards the right hand side
%
% requires the following functions from matlib:
% makeSubplotScalesEqual(), nancat(), conditionalPlot(), sq(),
% 
% nIter is the number of iterations across the parameter sweep (default = 100)
%
% it will %save the output in ./Data

%% simulate model

numTrials       = 100; % number of trials
itemsPerTrial   = ones(1,numTrials)*3; % number of items per trial
model           = SwapModel2D(); % 2D misbinding model
if nargin == 0
    nIter = 100;
end
nSteps = 6;
bias = linspace(-.2,.2,nSteps); % radial bias size

%% grid search of params
nParSteps = 11;
cols = [3:2:10];
g = linspace(0.01,.98,nParSteps);
b = linspace(0.01,.98,nParSteps);
sd = linspace(0.1,100,nParSteps);

params = [ kron(g, ones(1,nParSteps^2)) ;...
    kron(ones(1,nParSteps), kron(b, ones(1,nParSteps)));...
    kron(ones(1,nParSteps^2), sd)]';

params(sum(params(:,1:2),2) >1,:) = [];

params = repmat(params,nIter,1);

nParSets = length(params);
dimensions = [1366 768];
%%
[fitPars1, fitPars2,fitPars3] = deal(NaN(nParSets,3,nSteps));
fitPars3(:,4,:) = NaN;
targs = cell(nParSets,nSteps);
biasCoords = [dimensions(1); NaN]; % right edge of screen

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

            resps = simData.errors + simData.targets; % response coords
            diffs = resps - biasCoords; % distance to screen edge on each response
            diffs(isnan(diffs)) = 0; % set NaN to zero for later addition
            
            move = diffs .* bias(j); % size of distance to move
            
            resp2 = resps - move; % move the responses
            for iDim = 1:2 % replace any responses outside area
                resp2(resp2(:,iDim) > dimensions(iDim),iDim) = dimensions(iDim);
                resp2(resp2(:,iDim) <0,iDim) = 0;
            end
            
            err2 = resp2 - simData.targets; % apply bias
            
            simData2 = simData;
            simData2.errors = err2;
            simData2.biasCoords = biasCoords;
            fitPars2(i,:,j) = MLE(simData2, model);
            fitPars3(i,:,j) = MLE(simData2,WithRadialBias2D(model));
        end
        
    
end
t = toc
save('Data/MemToolbox2DSimPropBias.mat')
%%

simPars = params;

% reorder pars and calcualte a=1-g-B;
fitParsOrdered = [fitPars2(:,3,:), 1 - fitPars2(:,1,:) - fitPars2(:,2,:), fitPars2(:,2,:), fitPars2(:,1,:)];
simParsOrdered = repmat([simPars(:,3), 1 - simPars(:,1) - simPars(:,2), simPars(:,2), simPars(:,1)], 1,1,nSteps);
fitPars3Ordered = [fitPars3(:,3,:), 1 - fitPars3(:,1,:) - fitPars3(:,2,:), fitPars3(:,2,:), fitPars3(:,1,:), fitPars3(:,4,:)];

%%
figure()

m = ['ko';'bx';'g+';'r^'];
modelNames = {'1D','2D'};
parNames = {'\sigma','\alpha','\beta','\gamma'}; % standard deviation, target, misbind, guess
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

absMeanErr= abs(nancat(4, fitParsOrdered,fitPars3Ordered(:,1:4,:)) - simParsOrdered) ./ simParsOrdered;

figure(2);clf
c = get(gca,'ColorOrder');
labels = {'uncorrected','corrected'};
for i=1:2
    for k=1:2
        subplot(2,2,(k-1)*2+i)
        if k==1
            a = absMeanErr(:,:,:,i);
        else
            a = absMeanErr(5:605:nParSets,:,:,i);
        end

        h = errorBarPlot(permute(a(:,2:4,:),[1,3,2]),'area',1);
        hold on
        yyaxis right
        h2 = errorBarPlot(permute(a(:,1,:),[1,3,2]),'area',1);
        gc = gca;
        gc.YColor = c(1,:);
        h2{1,1}.Color = c(1,:);
        h2{1,1}.LineStyle = '-';
        h2{1,2}.FaceColor = c(1,:);
        for j = 1:3
            h{j,1}.Color = c(j+1,:);
            h{j,2}.FaceColor = c(j+1,:);
        end
        if i==2,    ylabel('SD abs error'), end
        yyaxis left
        gc = gca;
        gc.YColor = [0 0 0];


        xlim([1 nSteps])
        set(gca,'XTick',1:nSteps,'XTickLabel',bias)
        if i==1, ylabel('abs mean error'), end
        xlabel('radial bias')
        title(labels{i})
        box off
    end
end
legend([[h2{:,1}],[h{:,1}]],parNames,'Location','Best')
makeSubplotScalesEqual(2,2,1:2,[2 2], 0)
makeSubplotScalesEqual(2,2,3:4,[2 2], 0)
% saveas(figure(2), 'Figs/MemToolbox2DSimPropBias_1.jpg');
%%
figure(3)
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
SuperTitle('uncorrected')
%saveas(figure(3), 'Figs/MemToolbox2DSimPropBias_2.jpg');
%%
figure(4)
for j = 1:nSteps
    for i=1:4
        subplot(nSteps,4,(j-1)*4+i)
        s = simParsOrdered(:,i);
        
        f = fitPars3Ordered(:,i,j);
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
%saveas(figure(4), 'Figs/MemToolbox2DSimPropBias_3.jpg');

%% bias recovery

figure(5);clf;
errorBarPlot(sq(fitPars3Ordered(:,5,:)),'area',1,'plotargs',{'Marker','x'})
set(gca,'XTick',1:6, 'XTickLabel',bias, 'YTick',bias)
% axis([1 6 -.35 .35])
box off
xlabel('bias')
ylabel('recovered bias')
%saveas(figure(5), 'Figs/MemToolbox2DSimPropBias_4.jpg');

%% plot diff in pars vs pars
d = nancat(5, fitParsOrdered(:,1:4,:) - simParsOrdered, fitPars3Ordered(:,1:4,:) - simParsOrdered);
parNames = {'\sigma','\alpha','\beta','\gamma'};


figure(6);clf;
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
           ylabel(sprintf('%s\nrecovery error',labels{j}),'FontWeight','bold')
       end
    end
end
subplot(2,4,5), colorbar(gca,'North','Ticks',[0 .5 1],'TickLabels',[bias(1) 0 bias(end)])
makeSubplotScalesEqual(2,4,[2:4 6:8])
makeSubplotScalesEqual(2,4,[1 5])
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
xlabel('simulated parameter','FontWeight','normal')
%saveas(figure(6), 'Figs/MemToolbox2DSimPropBias_5.jpg');



%% bland altman

figure(7);clf
meanVals = nancat(4, fitParsOrdered(:,1:4,:) + simParsOrdered, fitPars3Ordered(:,1:4,:) + simParsOrdered) ./ 2;

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
           ylabel(sprintf('%s\nfit - sim',labels{j}));
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
subplot(2,4,8), colorbar(gca,'South','Ticks',[0 1], 'TickLabels',arrayfun(@num2str, bias([1 end]),'UniformOutput',0))
    
%saveas(figure(7), 'Figs/MemToolbox2DSimPropBias_6.jpg');

%% plot diff in pars vs pars
d = abs(nancat(5, fitParsOrdered(:,1:4,:) - simParsOrdered, fitPars3Ordered(:,1:4,:) - simParsOrdered));
parNames = {'\sigma','\alpha','\beta','\gamma'};


figure(8);clf;
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
xlabel('simulated parameter','FontWeight','normal')
%saveas(figure(8), 'Figs/MemToolbox2DSimPropBias_7.jpg');


%%

save('Data/MemToolbox2DSimPropBias.mat')