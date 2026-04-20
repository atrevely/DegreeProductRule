function [critJump, jumpStd, locMaxJump, fracFullyConn, maxClustMeanFld, maxClustSelfCrit] = PercMaxJumpMeans(numNodes, numChoices, sPointF, ePointF, sPointB, ePointB, numRuns, len, meancrit, cType, alpha)
% For a specified growth process (DPR or Achlioptas), performs a given
% number of growth realizations, parallelized, then averages together the
% evolution of the order parameter (largest cluster scaled by system size).
% Calculates the mean and standard deviation of the location of the largest
% jump in the order parameter.

% numNodes = number of nodes in the network.
% numChoices = number of edges competing at each timestep.
% sPointF = timestep where we'd like to start calculating the order
% parameter for regular percolation (forwards).
% ePointF = timestep where we're done calculating the order parameter for 
% regular percolation (forwards).
% sPointB = timestep where we'd like to start calculating the order
% parameter for reverse percolation.
% ePointB = timestep where we're done calculating the order parameter for 
% reverse percolation.
% numRuns = number of realizations to average together.
% len = total number of timesteps (before reversing the process).
% meancrit = critical point where we'd like to calculate maximum cluster.
% cType = type of network growth.
% alpha = paramete for interpolating between Erdos-Renyi and DPR.

% critJump = mean maximum jump in the order parameter.
% jumpStd = standard deviation of the maximum jump distribution.
% locMaxJump = timestep where the maximum jump occured for each
% realization.
% fracFullyConn = number of runs that reach full connectivity before
% terminating.
% maxClustMeanFld = size of the largest cluster at the mean field
% approximation of the critical point.
% maxClustSelfCrit = size of the largest cluster at the point in each
% individual run where the jump in the order parameter was largest.


% Initialize the max jump variables.
cJumpMaxF = zeros(1,numRuns);
cJumpMaxB = zeros(1,numRuns);
locMaxJumpF = zeros(1,numRuns);
locMaxJumpB = zeros(1,numRuns);

% Run realizations, storing the relevant numbers from each run.
passedSeed = randi(2^30,[1,numRuns]); % These need to be created OUTSIDE the parfor loop otherwise you'll get repeated seeds (thanks parfor).
parfor a = 1:numRuns
    rng(passedSeed(a));
    
    % Forward percolation. The 1e10 (and ~) is hard-coded because we don't really
    % want the distributions here (that's what the other code is for).
    if strcmp(cType,'PR') == 1
        [rEdge, cEdge, ~] = AchlioptasGrowthProcess(numNodes,numChoices,len,alpha,1e10);
    else
        [rEdge, cEdge, ~] = DPRGrowthProcess(numNodes,numChoices,len,alpha,1e10);
    end
    % Calculate order parameter forward, then store the largest jump.
    % Also store the size of the largest cluster at the largest jump.
    [fracConnF] = ClusterSizesRange(rEdge,cEdge,numNodes,sPointF,ePointF);
    [cJumpMaxF(a), locMaxJumpF(a)] = max(diff(fracConnF));
    maxClustMeanFld(a) = fracConnF(round(meancrit*numNodes) + 1 - sPointF);
    maxClustSelfCrit(a) = fracConnF(locMaxJumpF(a)+1);
    fConEnd(a) = fracConnF(length(fracConnF));
    
    % Reverse percolation.
    if strcmp(cType,'PR') == 1
        [rEdgeB, cEdgeB] = RExpPercAchlioptas(numNodes,numChoices,rEdge,cEdge);
    else
        [rEdgeB, cEdgeB] = RExpPercFast(numNodes,numChoices,rEdge,cEdge);
    end
    % Calculate order parameter backward, then store largest drop.
    [fracConnB] = ClusterSizesRange(rEdgeB,cEdgeB,numNodes,sPointB,ePointB);
    [cJumpMaxB(a), locMaxJumpB(a)] = max(diff(fracConnB));
end

% Add back in the steps that were skipped by windowing.
locMaxJumpF = locMaxJumpF + sPointF;
locMaxJumpB = locMaxJumpB + sPointB;

% Take means of forwards and backwards max jumps/drops.
critJump(1) = mean(cJumpMaxF);
critJump(2) = mean(cJumpMaxB);
jumpStd(1) = std(cJumpMaxF);
jumpStd(2) = std(cJumpMaxB);
locMaxJump(1,:) = locMaxJumpF;
locMaxJump(2,:) = locMaxJumpB;
% Return the fraction of runs that don't reach full connectedness.
fracFullyConn = nnz(fConEnd == 1)/numRuns;
end