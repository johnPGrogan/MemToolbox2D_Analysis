function o = WMTCalcBehVars(data)
% o = WMTCalcBehVars(data)
% calculates behavioural metrics for spatial WM tasks on simulated data
% e.g. mean target distance, nearest neighbour distance etc
% data is a structure with fields: distractors, errors, resps, targets, items, dimensions
% o is a structure containing all metrics

%%
distractors = data.distractors;% ./ pixelsPerCm(i);
errors = data.errors;

if ~isfield(data,'resps')
    if isfield(data, 'targets')
        data.resps = data.errors + data.targets;
    elseif isfield(data,'items') && isfield(data,'whichIsTestItem')
        iTest = data.whichIsTestItem;
        for i = 1:size(data.errors,2)
            data.resps(:,i) = data.items(iTest(i)*2-1:iTest(i)*2,i);
        end
        data.targets = data.resps - data.errors;
    else
        error('supply either data.resp, data.targets, data.items + data.whichIsTestItem');
    end
end

%%

targDists = sqrt (errors(1,:).^2 + errors(2,:).^2);
targDist = nanmean(targDists);

nItems = size(data.distractors,1)./2;
nnDists = NaN(nItems,size(data.errors,2));
for iDist = 1:(nItems)
    nnDists(iDist,:) = sqrt( (errors(1,:) - distractors((iDist*2-1),:)).^2 + (errors(2,:) - distractors(iDist*2,:)).^2 );
end
nnDists = min([targDists; nnDists]);
nnDist = nanmean(nnDists);
%%
nnMisb = nnDists < targDists; % trials where nn is nonTarget
[swapErr, swapErrMean] = deal(nnMisb');

annulusTol = 1.05 .* 41.269;%distance threshold for removing nn trials. 1.5vis deg in pixels
annulusTolMean = targDist; % mean target distance as threshold

swapErr(nnDists>annulusTol) = 0; % any outside thresh are non swaps
swapErrMean(nnDists>annulusTolMean) = 0; % any outside thresh are non swaps

stepSize = 1;%stepSize for circle
[circleX, circleY] = makeCircle(targDists', data.targets(1,:)', data.targets(2,:)', stepSize); % should this be centred on target?
xMin = 0;
xMax = data.dimensions(1);
yMin = 0;
yMax = data.dimensions(2);

circleX(circleX < xMin | circleX > xMax) = NaN;%remove any outside of boundaries
circleY(circleY < yMin | circleY > yMax) = NaN;

% find how many points are still within the screen
numPossibles = sum(~isnan(circleX) & ~isnan(circleY),2);

%see if each point is within annulusTol from any shape
withinShapeAnnuli = NaN(size(data.resps,2),361,nItems+1);
withinShapeAnnuliMean = withinShapeAnnuli;
for k = 1:nItems+1
    withinShapeAnnuli(:,:,k) = ((circleX - repmat(data.items(k*2-1,:)',1,361)).^2 + (circleY - repmat(data.items(k*2,:)',1,361)).^2) < annulusTol.^2;
    withinShapeAnnuliMean(:,:,k) = ((circleX - repmat(data.items(k*2-1,:)',1,361)).^2 + (circleY - repmat(data.items(k*2,:)',1,361)).^2) < annulusTolMean.^2;
end

totalWithins = sum(any(withinShapeAnnuli(:,:,2:end),3),2);%number of points within any shape's annuli
totalWithins(swapErr==0) = 0;%only keep this for trials with swap errors

probResp = totalWithins./numPossibles;
swapErrCorrected = swapErr - probResp;

swapErrors = nanmean(swapErr).*100;
swapErrorsCorrected = nanmean(swapErrCorrected).*100;

%mean threshold
totalWithinsMean = sum(any(withinShapeAnnuliMean,3),2);%number of points within any shape's annuli
totalWithinsMean(swapErrMean==0) = 0;%only keep this for trials with swap errors

probRespMean = totalWithinsMean./numPossibles;
swapErrMeanCorrected = swapErrMean - probRespMean;

swapErrorsMean = nanmean(swapErrMean).*100;
swapErrorsMeanCorrected = nanmean(swapErrMeanCorrected).*100;


o = workspace2struct();
end

function [circleX,circleY] = makeCircle(radius,centreX,centreY,stepSize)
% draws a circle around a point
% [x,y] = makeCircle(radius,centreX,centreY,stepSize)

    theta = -180:stepSize:180;%angle of cirlce in stepSizes
    thetaRadians = theta * (pi/180);%convert to radians

    %now calc circle coordinates, centered on the responseLocation
    circleX = radius*sin(thetaRadians) + repmat(centreX,1,length(thetaRadians));
    circleY = radius*cos(thetaRadians) + repmat(centreY,1,length(thetaRadians));
    
end