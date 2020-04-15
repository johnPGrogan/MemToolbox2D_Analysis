function MemToolbox2DSimResampling(nIter)
% MemToolbox2DSimResampling(nIter)
% simulation with edge-resampling and edge-constraining methods for responses outside the screen
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
models          = {SwapModel2D(), SwapModel2D(1)}; % 1D and 2D misbinding models
nPPs            = 100; % number of pps to simulate
if nargin == 0
    nIter = 100;
end
nSteps = 11;
g = linspace(0.01,.98,nSteps);
b = linspace(0.01,.98,nSteps);
sd = linspace(0.1,100,nSteps);

params = [ kron(g, ones(1,nSteps^2)) ;...
    kron(ones(1,nSteps), kron(b, ones(1,nSteps)));...
    kron(ones(1,nSteps^2), sd)]';
params(sum(params(:,1:2),2) >1,:) = [];

params = repmat(params,nIter,1);

nParSets = size(params, 1);
%% fit model

fitParsBound = cell(nParSets,1); %preallocate
llBound      = cell(nParSets,1); % log likelihoods
fitParsResamp = fitParsBound; %preallocate
llResamp      = llBound; % log likelihoods

tic
parfor i = 1:nParSets
    disp(i)
        

    % get locations
    simulatedData2D{i} = GenerateDisplays2D(numTrials, itemsPerTrial); 
    
    simulatedData2DBound{i} = simulatedData2D{i}; % make a copy for the boundary simulation
    
    % simulate responses - with resampling
     simulatedData2D{i}.errors = SampleFromModel2D(models{1}, params(i,:),...
        [1 numTrials], simulatedData2D{i}); % simulate responses
     
     % fit swap model with maximum likelihood
    [fitParsResamp{i,1}, llResamp{i,1}] = MLE(simulatedData2D{i}, models{1});
   

    % with boundary

%     simulatedData2DBound{i}.boundary = 1;
    simulatedData2DBound{i}.errors = SampleFromModel2D(models{2}, params(i,:),...
        [1 numTrials], simulatedData2DBound{i}); % simulate responses
    
    % fit swap model with maximum likelihood
    [fitParsBound{i,1}, llBound{i,1}] = MLE(simulatedData2DBound{i}, models{2});

    
end
toc
save('Data/MemToolbox2DSimResampling.mat')
%% plot params from swap model

fitPars = cat(3, cat(1,fitParsBound{:}), cat(1,fitParsResamp{:}));

% reorder pars and calcualte a=1-g-B;
fitParsOrdered = [fitPars(:,3,:), 1 - fitPars(:,1,:) - fitPars(:,2,:), fitPars(:,2,:), fitPars(:,1,:)];
simParsOrdered = repmat([params(:,3), 1 - params(:,1) - params(:,2), params(:,2), params(:,1)], 1,1,2);

%%
modelNames = {'resampling','constraining'};
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
            axis([0 max(sd) 0 max(sd)])
        else
            axis([0 1 0 1])
        end
    end
end
%saveas(figure(1), 'Figs/MemToolbox2DSimResampling_1.jpg');
%% heatmaps
figure()

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
figure(3)
for j = 1:2
    for i=1:4
        subplot(2,4,(j-1)*4+i)
        s = simParsOrdered(:,i,j);

        
        f = fitParsOrdered(:,i,j);
        fp = groupMeans(f,1,s,'dim');
        quantPlot(unique(s),fp);
        hold on
        line([0 max(sd)],[0 max(sd)],'Color','k','LineStyle','--')
        if j==1
            title(parNames{i})
        else
            xlabel('sim pars')
        end
        if i==1
            ylabel(modelNames{j})
            axis([0 max(sd) 0 max(sd)])
        else
            axis([0 1 0 1])
        end
    end
end
%saveas(figure(3), 'Figs/MemToolbox2DSimResampling_2.jpg');

%% plot diff in pars vs pars
d = cat(4, fitParsOrdered - simParsOrdered, abs(fitParsOrdered - simParsOrdered));
for i = 1:4
    [~,diffP(i)] = ttest(d(:,i,1), d(:,i,2));
end
labels = {'Recovery error','Abs( recovery error )'};
modelNames = {'Resampling','Constraining'};
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
       if j==1, title(parNames{i}); end
       if j==2,    xlabel(parNames{i},'FontWeight','bold'), end
       if i==1, ylabel(labels{j},'FontWeight','bold'), end
       box off
       x = unique(round(simParsOrdered1(:,i,1), 0+(i>1)));
       set(gca,'XTick',1:5:11,'XTickLabel',x(1:5:end))
       xlim([1 11])
       if i==1; ylim([-20 20]); set(gca,'YTick',linspace(ylim(1), ylim(2), 3));
       else; ylim([-.06 .06]); set(gca,'YTick', linspace(ylim(1), ylim(2), 3));end
    end
end
makeSubplotScalesEqual(2,4,[2:4])
makeSubplotScalesEqual(2,4,[6:8])
legend(h(:,1),modelNames,'Location','Best')
h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('Simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;
%saveas(figure(4), 'Figs/MemToolbox2DSimResampling_3.jpg');
%%
save('Data/MemToolbox2DSimResampling.mat')