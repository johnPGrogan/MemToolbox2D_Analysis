function MemToolbox2DBehMetrics
% MemToolbox2DBehMetrics

% clear
% close all
load('./Data/MemToolbox2DSim.mat', 'params','simulatedData2D','fitPars2D','simParsOrdered');

%% remove invalid pars

toRemove = sum(params(:,1:2),2) >= 1;
simulatedData2D(toRemove) = [];

%% calculate beh metrics for simulated data

n = size(simulatedData2D,1);

behVarNames = {'targDist', 'nnDist', 'swapErrors','swapErrorsCorrected',...
    'swapErrorsMean', 'swapErrorsMeanCorrected'};

parfor i = 1:n
    if mod(i,100)==0; disp(i); end
    
    behVars = WMTCalcBehVars(simulatedData2D{i}); % calculate metrics
    
	% extract
    targDist(i,1) = behVars.targDist;
    nnDist(i,1) = behVars.nnDist;
    swapErrors(i,1) = behVars.swapErrors;
    swapErrorsCorr(i,1) = behVars.swapErrorsCorrected;
    swapErrorsMean(i,1) = behVars.swapErrorsMean;
    swapErrorsMeanCorr(i,1) = behVars.swapErrorsMeanCorrected;
end

clear simulatedData2D % clear to reduce memory load

%% put into mat

behVarMat = [targDist, nnDist, targDist-nnDist, swapErrors, swapErrorsCorr, swapErrorsMean, swapErrorsMeanCorr];


%% only use 10 iterations for memory reasons
nIter = min(size(fitPars2D,3), 10); % max 10 iterations of par sweep

fitPars2D(toRemove,:,:) = [];
simParsOrdered(toRemove,:) = [];
simParsOrdered = repmat(simParsOrdered, 1, 1, nIter);

fitPars = fitPars2D;
fitParsOrdered = [fitPars(:,3,:), 1 - fitPars(:,1,:) - fitPars(:,2,:), fitPars(:,2,:), fitPars(:,1,:)];

fitParsOrdered = fitParsOrdered(:,:,1:nIter);
behVarMat = repmat(behVarMat, 1, 1, nIter);



%% plot beh vars against sim pars
parNames = {'\sigma','\alpha','\beta','\gamma'};
varNamesShort = {'TD','NN','TD-NN','SE','SEC','SEm','SECm'};
varNames = {'targDist', 'nearestNeighbour','targDist - nearestNeighbour','swapError','swapErrorCorr','swapErrorMeanThresh','swapErrorMeanThreshCorr'};

figure()
for iVar = 1:7
    for iPar = 1:4
        subplot(7,4,(iVar-1)*4+iPar)
        plot(col(simParsOrdered(:,iPar,:)),col(behVarMat(:,iVar,:)),'.');
        [correls(iVar,iPar,1),correls(iVar,iPar,2)] = corr(col(simParsOrdered(:,iPar,:)),...
            col(behVarMat(:,iVar,:)),'type','pearson');
        if correls(iVar,iPar,2) < .05
            l = lsline;
            l.Color = [1 0 0];
            l.LineWidth = 2;
        end
        if iPar==1
            ylabel(varNames{iVar})
        else
            xlim([0 1])
        end
        if iVar==1
            title(parNames{iPar})
        elseif iVar==7
            xlabel(parNames{iPar})
        end
        box off
    end
end

%%

figure()
for iVar = 1:7
    for iPar = 1:4
        subplot(4,7,(iPar-1)*7+iVar)
        plot(col(behVarMat(:,iVar,:)), col(fitParsOrdered(:,iPar,:)),'.');
        [r,p] = corrcoef(col(behVarMat(:,iVar,:)), col(fitParsOrdered(:,iPar,:)));
        correls2(iPar, iVar, 1) = r(2);
        correls2(iPar, iVar, 2) = p(2);
        if p(2) < .05
            l = lsline;
            l.Color = [1 0 0];
            l.LineWidth = 2;
        end
        if iPar==1
            %             title(varNames{iVar})
        else
            ylim([0 1])
            if iPar==4
                xlabel(varNames{iVar})
            end
        end
        if iVar==1
            ylabel(parNames{iPar})
        end
        box off
    end
end

%%

figure()
for iVar = 1:4
    for iPar = 1:4
        subplot(4,4,(iPar-1)*4+iVar)
        plot(col(simParsOrdered(:,iPar,:)), col(fitParsOrdered(:,iVar,:)),'.');
        [r,p] = corrcoef(col(simParsOrdered(:,iPar,:)), col(fitParsOrdered(:,iVar,:)));
        correls3(iPar, iVar, 1) = r(2);
        correls3(iPar, iVar, 2) = p(2);
        if p(2) < .05
            l = lsline;
            l.Color = [1 0 0];
            l.LineWidth = 2;
        end
        if iPar==1
            %             title(par{iVar})
        elseif iPar==4
            xlabel(parNames{iVar})
        end
        if iVar==1
            ylabel(parNames{iPar})
        end
        box off
    end
end

%%

figure();
inds = [1 2; 1 4; 2 1; 3 3; 4 2];
correls4 = NaN(2,5,2);
for i = 1:5
    subplot(2,5,i)
    plot( col(simParsOrdered(:,inds(i,1),:)),col(behVarMat(:,inds(i,2),:)),'.');
    [correls4(1,i,1),correls4(1,i,2)] = corr(col(simParsOrdered(:,inds(i,1),:)), ...
        col(behVarMat(:,inds(i,2),:)),'type','pearson');
    if p(2) < .05
        l = lsline;
        l.Color = [1 0 0];
        l.LineWidth = 2;
    end
    if i==1
        ylabel(sprintf('%s',varNames{inds(i,2)}))
    else
        ylabel(varNames{inds(i,2) })
    end
    if i>2
        xlim([0 1])
    end
    xlabel(parNames{inds(i,1)},'FontWeight','bold')
    
    box off
end

for iPar = 1:4

    subplot(2,5,iPar+6)
    
    plot(col(simParsOrdered(:,iPar,:)), col(fitParsOrdered(:,iPar,:)),'.');
    [correls4(2,iPar+1,1),correls4(2,iPar+1,2)] = corr(col(simParsOrdered(:,iPar,:)), ...
        col(fitParsOrdered(:,iPar,:)),'type','Spearman');

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
x = xlabel('simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;


%% quant plot

percs = [0 0.2 2.3 15.9 50 84.1 97.7 99.8 100];
figure(5);clf

for i = 1:5
    subplot(2,5,i)
    s = round(simParsOrdered(:,inds(i,1)),1);
    f = behVarMat(:,inds(i,2));
    
    fp = groupMeans(f,1,s,'dim');
    quantPlot(unique(s),fp,percs);
    
    if i==1
        ylabel(sprintf('%s',varNames{inds(i,2)}))
    else
        ylabel(varNames{inds(i,2) })
    end
    if i>2
        xlim([0 1])
    end
    xlabel(parNames{inds(i,1)},'FontWeight','bold')
    
    box off
end

for iPar = 1:4
    subplot(2,5,iPar+6)
    s = round(simParsOrdered(:,iPar),1);
    f = fitParsOrdered(:,iPar);
    
    fp = groupMeans(f,1,s,'dim');
    quantPlot(unique(s),fp,percs);
        
    if iPar==1
        ylim([0 100])
        ylabel(parNames{iPar},'FontWeight','bold')
    else
        ylabel(parNames{iPar},'FontWeight','bold')
        xlim([0 1])
    end
    xlabel(parNames{iPar},'FontWeight','bold')
    

    box off
end

h = axes('visible','off'); % super X and Y labels
h.XLabel.Visible = 'on';
x = xlabel('simulated parameter','FontWeight','normal');
x.Position(2) = x.Position(2)*1.2;


%%

figure(6);clf
for iVar = 1:7
    for iPar = 1:4
        subplot(7,4,(iVar-1)*4+iPar)
        s = round(simParsOrdered(:,iPar),1);
        f = behVarMat(:,iVar);

        fp = groupMeans(f,1,s,'dim');
        quantPlot(unique(s),fp,percs);

        if iPar==1
            ylabel(sprintf('%s',varNamesShort{iVar}))

        end
        if iPar>2
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

%%

save('./Data/MemToolbox2DBehMetrics.mat','simParsOrdered','fitParsOrdered','behVarMat','p')
end