function MemToolbox2DSim(nIter)
% MemToolbox2DSim(nIter)
% simulate and fit 2D swap model over a large range of parameters
% plot figures looking at quality of fit
%
% requires the following functions from matlib:
% makeSubplotScalesEqual(), nancat(), conditionalPlot(), sq(),
% 
% nIter is the number of iterations across the parameter sweep (default = 100)
%
% it will save the output in ./Data
%

%%

numTrials       = 100; % number of trials
itemsPerTrial   = ones(1,numTrials)*3; % number of items per trial
model           = SwapModel2D(); % 2D misbinding model
nPPs            = 100; % number of pps to simulate

if nargin==0
    nIter = 2;
end
%% grid search of params
nSteps = 11;
cols = [3:2:10];
g = linspace(0.0,1,nSteps);
b = linspace(0,1,nSteps);
sd = linspace(0.01,100,nSteps);


params = [ kron(g, ones(1,nSteps^2)) ;...
    kron(ones(1,nSteps), kron(b, ones(1,nSteps)));...
    kron(ones(1,nSteps^2), sd)]';

nParSets = length(params);
%%
fitPars2D = NaN(nParSets,3,nIter);
ll = NaN(nParSets,1,nIter);
simulatedData2D = cell(nParSets,1); 
tic
parfor i=1:nParSets
    disp(i)
    
    if sum(params(i,1:2)) <= 1 % can't simulate if 
        for j=1:nIter
            % get locations
            simData = GenerateDisplays2D(numTrials, itemsPerTrial);
            
            % simulate responses
            simData.errors = SampleFromModel2D(model, params(i,:),...
                [1 numTrials], simData); % simulate responses
            
            
            % fit swap model with maximum likelihood
            [fitPars2D(i,:,j),ll(i,1,j)] = MLE(simData, model);
            
            simulatedData2D{i} = simData;
        end
    end
end
toc
save('Data/MemToolbox2DSim.mat')
%%

params2 = permute(reshape(repmat(params,1,1,nIter),nSteps,nSteps,nSteps,3,nIter),[3,4,2,1,5]);
fitPars = permute(reshape(fitPars2D,nSteps,nSteps,nSteps,3,nIter),[3,4,2,1,5]);
nll = - permute(reshape(ll,nSteps,nSteps,nSteps,1,nIter),[3,4,2,1,5]);

fitParsOrdered = [fitPars2D(:,3,:), 1 - fitPars2D(:,1,:) - fitPars2D(:,2,:), fitPars2D(:,2,:), fitPars2D(:,1,:)];
simParsOrdered = repmat([params(:,3), 1 - params(:,1) - params(:,2), params(:,2), params(:,1)], 1,1,1);

%% plot nll over par space
figure(1);clf
scatter3(col(params2(:,1,:,:,1)), col(params2(:,2,:,:,1)), col(params2(:,3,:,:,1)),50, col(nanmean(nll,5)),'.');
colormap('jet')
cb = colorbar;
xlabel('g')
ylabel('b')
zlabel('sd')
cb.Label.String = 'nll';
caxis([1200 1400])

%% mean absolute error across params

absMeanErr = squeeze(nanmean(abs(params2 - fitPars) ./ params2,[2,5]));

figure(2);clf
scatter3(col(params2(:,1,:,:,1)), col(params2(:,2,:,:,1)), col(params2(:,3,:,:,1)),50, col(absMeanErr),'.');
colormap('jet')
cb = colorbar;
xlabel('g')
ylabel('b')
zlabel('sd')
cb.Label.String = 'abs mean err';
caxis([0 1])

%% guess
legendLab = {'y=x','median','±1SD','±2SD','±3SD','±4SD'};
figure(3);clf
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
            elseif iB==n
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
xlabel('simulated \gamma value','FontWeight','bold')
h.YLabel.Visible = 'on';
y = ylabel('recovered \gamma value','FontWeight','bold');
y.Position(1) = y.Position(1) *1.2;

%saveas(figure(3), 'Figs/MemToolbox2DSim_1.jpg');
%% misbind
figure(4);clf
inds = [];
for iG = 1:n
    for iSD=1:n
        if any(~isnan(col(fitPars(cols(iG),2,:,cols(iSD),:))))
            subplot(n,n,(iG-1)*n+iSD)
            p = quantPlot(sq(params2(cols(iG), 2,:, cols(iSD),1)), sq(fitPars(cols(iG),2,:,cols(iSD),:)));
            hold on
            line([0 1],[0 1],'Color','k','LineStyle','--')
            if iG==1
                title(sprintf('\\sigma=%.1f',sd(cols(iSD))),'FontWeight','normal')
            elseif iG==n
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

%saveas(figure(4), 'Figs/MemToolbox2DSim_2.jpg');
%% SD
figure(5);clf
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
            elseif iG==n
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

%saveas(figure(5), 'Figs/MemToolbox2DSim_3.jpg');


%% plot diff in pars vs pars
f = load('./Data/MemToolbox1DSim.mat', 'fitParsOrdered','simParsOrdered');
fitPars1DOrdered = f.fitParsOrdered;
simPars1DOrdered = f.simParsOrdered;
d = nancat(4, fitPars1DOrdered - simPars1DOrdered, fitParsOrdered - simParsOrdered);

modelNames = {'1D','2D'};
parNames = {'\sigma','\alpha','\beta','\gamma'};

figure(6);clf;
simParsOrdered1(:,:,1) = round(simPars1DOrdered,2,'significant');
simParsOrdered1(:,:,2) = round(simParsOrdered,2,'significant');
simParsOrdered1(simParsOrdered1<0) = NaN;
c  = [1 0 0; 1 .4 0; 1 1 0;0 1 0; 0 1 1; 0 0 1;];
for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
%        set(gca,'ColorOrder',c);
       diff = reshape(groupMeans(d(:,i,:,j),1, repmat(simParsOrdered1(:,i,j),1,1,nIter),'dim'),11,[]);
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
xlabel('simulated parameter','FontWeight','normal')

%saveas(figure(6), 'Figs/MemToolbox2DSim_4.jpg');

%% bland altman

figure(7);clf
meanVals = nancat(4, fitPars1DOrdered + simPars1DOrdered, fitParsOrdered + simParsOrdered) ./ 2;
simParsOrdered1 = round(meanVals,2,'significant');

for j = 1:2
    for i = 1:4
       subplot(2,4,(j-1)*4+i)
        for k = 1
            conditionalPlot(sq(meanVals(:,i,:,j)), sq(d(:,i,:,j)));%,[],'color',c(k,:));
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
    
%saveas(figure(7), 'Figs/MemToolbox2DSim_5.jpg');


%%
save('Data/MemToolbox2DSim.mat')