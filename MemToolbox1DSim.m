function MemToolbox1DSim(nIter)
% MemToolbox1DSim(nIter)
% simulate and fit 1D swap model over a large range of parameters
% plot figures looking at quality of fit
%
% requires the following functions from matlib:
% makeSubplotScalesEqual(), nancat(), conditionalPlot(), sq(),
% 
% nIter is the number of iterations across the parameter sweep (default = 100)
%
% it will save the output in ./Data
%


close all

%% simulate and look at params

numTrials       = 100; % number of trials
itemsPerTrial   = ones(1,numTrials)*3; % number of items per trial
model           = SwapModel(); % 1D and 2D misbinding models
nPPs            = 100; % number of pps to simulate
if nargin == 0
    nIter = 100;
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
fitPars1 = NaN(nParSets,3,nIter);
ll = NaN(nParSets,1,nIter);
tic
parfor i=1:nParSets
    disp(i)
    
    if sum(params(i,1:2)) <= 1 % can't simulate if 
        for j=1:nIter

            % get locations
            simData = GenerateDisplays(numTrials, itemsPerTrial);

            % simulate responses
            simData.errors = SampleFromModel(model, params(i,:),...
                [1 numTrials], simData); % simulate responses

            % fit swap model with maximum likelihood
            [fitPars1(i,:,j),ll(i,1,j)] = MLE(simData, model);

        end       
    end
end
toc
% save('Data/MemToolbox1DSim.mat')
%%

params2 = permute(reshape(repmat(params,1,1,nIter),nSteps,nSteps,nSteps,3,nIter),[3,4,2,1,5]);
fitPars = permute(reshape(fitPars1,nSteps,nSteps,nSteps,3,nIter),[3,4,2,1,5]);
nll = - permute(reshape(ll,nSteps,nSteps,nSteps,1,nIter),[3,4,2,1,5]);

fitParsOrdered = [fitPars1(:,3,:), 1 - fitPars1(:,1,:) - fitPars1(:,2,:), fitPars1(:,2,:), fitPars1(:,1,:)];
simParsOrdered = repmat([params(:,3), 1 - params(:,1) - params(:,2), params(:,2), params(:,1)], 1,1,1);

%% plot nll over par space
figure(1)

scatter3(col(params2(:,1,:,:,1)), col(params2(:,2,:,:,1)), col(params2(:,3,:,:,1)),50, col(nanmean(nll,5)),'.');
colormap('jet')
cb = colorbar;
xlabel('g')
ylabel('b')
zlabel('sd')
cb.Label.String = 'nll';

%% mean absolute error across params

absMeanErr = squeeze(nanmean(abs(params2 - fitPars) ./ params2,[2,5]));

figure(2)

scatter3(col(params2(:,1,:,:,1)), col(params2(:,2,:,:,1)), col(params2(:,3,:,:,1)),50, col(absMeanErr),'.');
colormap('jet')
cb = colorbar;
xlabel('g')
ylabel('b')
zlabel('sd')
cb.Label.String = 'abs mean err';
caxis([0 1])  

%% take some slices across par space

%fix 2 pars, look at correl of other par
legendLab = {'y=x','median','±1SD','±2SD','±3SD','±4SD'};

% guess
figure(3);clf
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
        elseif iB==n
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
xlabel('simulated \gamma value','FontWeight','bold')
h.YLabel.Visible = 'on';
y = ylabel('recovered \gamma value','FontWeight','bold');
y.Position(1) = y.Position(1) *1.2;
%saveas(figure(3), 'Figs/MemToolbox1DSim_1.jpg');
%% misbind
figure(4);clf
for iG = 1:n
    for iSD=1:n
        subplot(n,n,(iG-1)*n+iSD)
        p = quantPlot(sq(params2(cols(iG), 2,:, cols(iSD),1)), sq(fitPars(cols(iG),2,:,cols(iSD),:)));
        hold on
        line([0 1], [0 1], 'Color','k','LineStyle','--')
        if iG==1
            title(sprintf('\\sigma=%.1f',sd(colsSD(iSD))), 'FontWeight','normal')
        elseif iG==n
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
xlabel('simulated \beta value','FontWeight','bold')
h.YLabel.Visible = 'on';
y = ylabel('recovered \beta value','FontWeight','bold');
y.Position(1) = y.Position(1) *1.2;

%saveas(figure(4), 'Figs/MemToolbox1DSim_2.jpg');
%% SD
figure(5);clf
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
xlabel('simulated \gamma value','FontWeight','bold')
h.YLabel.Visible = 'on';
y = ylabel('recovered \gamma value','FontWeight','bold');
y.Position(1) = y.Position(1) *1.2;
%saveas(figure(5), 'Figs/MemToolbox1DSim_3.jpg');
%%
save('Data/MemToolbox1DSim.mat')