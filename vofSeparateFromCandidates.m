function vofIndices = vofSeparateFromCandidates(fg,paramImg,mfsFile)
% vofSeparateFromCandidates separates the VOF from a candidate set of
% vertical streamlines (fg), based on the local information in a T1 map
% (paramImg) or other map. Specifically, it finds the point where the
% change in T1 median values is maximal.
% 
% INPUT
% =====
% * fg - a fiber group of vertical streamlines (output of ____)
% * paramImg - a nifti image of a quantitative T1 map (of any resolution),
%              aligned to the diffusion data
% * mfsFile - a nifti representation of the MFS label from FreeSurfer (see _____).
%
% OUTPUT
% ======
% * vofIndices - a list of VOF streamline indices from the input fg

%% Set slgorithm parameters
n = 15; % Divide the stremalines to bins of size (yMax-yMin)/n;
mergeBinsWithFewStreamlines = true; % Merge consecutive streamlines if they have too few streamlines
minimalStreamlinesNum = 35;
numNodes = 6; % Number of nodes to remove from streamline ends
if ~exist('mfsFile','var') || isempty(mfsFile)
    mfsFlag = false;
else
    mfsFlag = true;
end

%% Trim the first and last nodes of each streamline, to avoid partial volume effect
fg.fibers = cellfun(@(x) x(:,numNodes+1:end-numNodes), fg.fibers,'UniformOutput',false);

%% Extract the spatial position and the paramMdn (e.g., T1-Mdn) of each streamline
y = cellfun(@(x) median(x(2,:)), fg.fibers)';
perPointFlag = 0;
fg = dtiCreateQuenchStats(fg, 'param_median', 'param', perPointFlag, paramImg, 'nanmedian', 1); % Yhis is a version of Vistasoft's dtiCreateQuenchStats, where in line 218 a "nanstd" options exists
paramMdn = fgGetParams(fg,'param_median');

%% (optional) Ignore streamlines with too much extrapolated B1+
noB1prcnt = fgGetParams(fg,'noB1prcnt');
if ~all(noB1prcnt==1)
    B1PrcntThresh = 0.2; % For finding the separating border, ignore streamlines with more than 20% extrapolated B1+
    indices = find(noB1prcnt<=B1PrcntThresh);
    y = y(indices);
    paramMdn = paramMdn(indices);
end

%% Set the relevant bins along the posterior-anterior axis
binWidth = (max(y)-min(y))/n;
prctOvrlp = 0.75; % The percent overlap between consecutive bins
step = binWidth*(1-prctOvrlp); % 0.75 overlap is 0.25 step size
binEdges = [min(y), min(y)+binWidth];
bI = 2;
while true
    if binEdges(bI-1,1) + step + binWidth > max(y)
        break
    end
    binEdges(bI,1) = binEdges(bI-1,1) + step;
    binEdges(bI,2) = binEdges(bI,1) + binWidth;
    bI = bI + 1;
end

smallBins = [];
if mergeBinsWithFewStreamlines
    for bI = 1:size(binEdges,1)
        indices = find(y>=binEdges(bI,1) & y<=binEdges(bI,2));
        if length(indices)<35 && bI<size(binEdges,1)% Ignore regions with too few streamlines. 35-40 seems to work fine
            smallBins(bI) = 1;
            binEdges(bI+1,1) = binEdges(bI,1); % attach this bin to the next one
        end
    end
    binEdges(find(smallBins),:) = [];
end

% Exclude the first and last 2 bins. The border can't be there, and these
% regions are prone to edge artifact
ignoreFirstAndLastBins = 0;
if ignoreFirstAndLastBins
    binEdges([1,2,end-1:end],:) = [];
end

%% Find median T1 in each bin
paramMdnBins = nan(1,size(binEdges,1));
for ii = 1:size(binEdges,1)
    indices = find(y>=binEdges(ii,1) & y<=binEdges(ii,2));
    if length(indices)<35
        continue
    end
    paramMdnBins(ii) = nanmedian(paramMdn(indices));
end

% Fix cases of NaN in some of the bins
notnanIndices = find(~isnan(paramMdnBins));
tmpParamMdnBins = paramMdnBins(notnanIndices);
paramMdnDiff = diff(tmpParamMdnBins);
[~,maxDiffIndex] = sort(paramMdnDiff,'descend');
if isempty(maxDiffIndex) | paramMdnDiff(maxDiffIndex)<0
    vofIndices = 1:length(fg.fibers);
    warning('Could not find any sharp increae in T1 value. Returning all streamlines');
    return
end

if mfsFlag
    % Find the posterior end of MFS (mfsY)
    mfs = readFileNifti(mfsFile);
    mfsIndices = find(mfs.data);
    [~,yMfsVox,~] = ind2sub(size(mfs.data), mfsIndices); % y vales of MFS in terms of voxels
    yMfsLims = minmax(yMfsVox'); % y values of MFS in terms of voxels
    tmp = mrAnatXformCoords(mfs.qto_xyz, [0,yMfsLims(1),0; 0, yMfsLims(2), 0]);
    mfsY = tmp(1,2); % y values of MFS in terms of image coordiantes
    
    % Choose between the two sharpest increases in T1, such that the VOF
    % terminates closest to the MFS
    indices_for_borders{1} = find(y<=binEdges(maxDiffIndex(1)+1,1)); % Find the indices of streamlines posterior to the 1st candidate-border bin
    indices_for_borders{2} = find(y<=binEdges(maxDiffIndex(2)+1,1)); % Find the indices of streamlines posterior to the 2nd candidate-border bin
    
    y1Tmp = y(indices_for_borders{1});
    y2Tmp = y(indices_for_borders{2});
    yPrctile(1) = prctile(y1Tmp,100); % One could play with the percentile here, and see if this leads to a more robust result
    yPrctile(2) = prctile(y2Tmp,100);
    mfsVofDist(1) = abs(yPrctile(1)-mfsY);
    mfsVofDist(2) = abs(yPrctile(2)-mfsY);
    
    [~,closeToMfsIdx] = min(mfsVofDist);
    maxDiffIndex = maxDiffIndex(closeToMfsIdx);
else
    maxDiffIndex = maxDiffIndex(1);
end

%% Return the VOF indices
maxDiffIndex = floor((notnanIndices(maxDiffIndex)+notnanIndices(maxDiffIndex+1))/2); % In case there are many NaNs after the border point, go to the middle of them
yThreshold = binEdges(maxDiffIndex+1,1);
y = cellfun(@(x) median(x(2,:)), fg.fibers)'; % Extract y again, in case some of the streamlines were removed above
vofIndices = find(y<=yThreshold);

