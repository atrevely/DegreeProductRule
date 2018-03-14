function [fConFMeans, clusterHistsComb] = PercAveragesRange(numNodes, numChoices, sPoint, ePoint, numRuns, cType, alpha, distPts)
% Runs many parallelized DPR Growth Processes or Achlioptas Processes, then averages the
% results and returns a vector of average size of the largest cluster at
% each timestep, as well as cluster size distributions at the specified
% points.

% numNodes = number of nodes in the network.
% numChoices = number of edges competing at each timestep.
% sPoint = timestep where we'd like to start calculating the order
% parameter.
% ePoint = timestep where we're done calculating the order parameter (and
% run stops).
% numRuns = number of runs we'd like to average together.
% cType = type of network growth.
% alpha = paramete for interpolating between Erdos-Renyi and DPR.
% distPts = points for calculating the cluster size distribution.

% fConFMeans = average of the order parameter at each timestep during the
% runs.
% clusterHistsComb = combined cluster size distribution at specified
% times across all the runs.

% Initialize the means and histogram cell.
fConF = [];
clusterHists = {};

% Run realizations, concatenating runs one after the other.
passedSeed = randi(2^30,[1,numRuns]); % These need to be created OUTSIDE the parfor loop otherwise you'll get repeated seeds (thanks parfor).
parfor a = 1:numRuns
    rng(passedSeed(a));
    % Forward percolation.
    if strcmp(cType,'PR') == 1
        [rEdge, cEdge] = AchlioptasGrowthProcess(numNodes,numChoices,ePoint,alpha); % If we're doing Achlioptas Process
    else
        [rEdge, cEdge] = DPRGrowthProcess(numNodes,numChoices,ePoint,alpha); % If we're doing DPR Process.
    end
    % Calculate order parameter and histograms.
    [fracConnF, cHists] = clstSizeRng(rEdge,cEdge,numNodes,sPoint,ePoint,distPts);
    fracConnF = fracConnF';
    cHists = cHists';
    fConF = [fConF, fracConnF];
    clusterHists = [clusterHists, cHists];
end

% Average all of the forward runs.
fConFMeans = sum(fConF,2)/numRuns;

% Pad cluster size histograms to make them equal length across runs, then combine them.
for a = 1:size(clusterHists,1)
    hstLens = cellfun(@length, clusterHists(a,:));
    for b = 1:length(clusterHists(a,:))
        clusterHists{a,b} = padarray(clusterHists{a,b},[0, max(hstLens) - hstLens(b)],'post');
    end
    clusterHistsComb{a} = sum(cat(1,clusterHists{a,:}));
end
end