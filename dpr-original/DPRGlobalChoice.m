function [maxClust, maxJump, nofinish] = DPRGlobalChoice(numNodes)
% Performs network growth using the DPR Process with global choice, meaning
% that every edge is considered for addition at each timestep. This is a
% bit of a tricky process, and is streamlined by recognizing some
% symmetries and simplifications in how the process eveolves.

% numNodes = number of nodes in the network.

% maxClust = size of the largest cluster at each timestep.
% maxJump = maximum jump in the size of the largest cluster at any single
% timestep.
% nofinish = indicates a run where the process did not finish connecting
% the entire network. These are exceedingly rare but require additional
% attention.

% Initialize cluster and loop variables.
clusters = 2*ones(1,numNodes/2);
loopNum = 1;
m = 1;
maxClust = zeros(1,numNodes/2);
maxClust(1) = 2/numNodes; % We know that we start with a cluster of size 2.
closedLoops = [];
nofinish = 0;
rng('shuffle');

% Here, nodes will connect and snake, eventually forming loops.
for q = 1:(numNodes/2-1)
    
    % Pull two random clusters to merge.
    clustChoices = sort(randperm(numNodes/2 - q + 1,2));
    
    % Figure out if it's a self-loop. Otherwise proceed as normal.
    if (clusters(clustChoices(1)) ~= 2) && (randi(numNodes - 2*q + 1) == 1)
        
        % This is the procedure for a closing a snake into loop.
        closedLoops(loopNum) = clusters(clustChoices(1));
        clusters([clustChoices(1) numNodes/2-q+1]) = clusters([numNodes/2-q+1 clustChoices(1)]);
        loopNum = loopNum + 1;
        maxClust(q + 1) = maxClust(q);
        
    else
        
        % Sum the two clusters in one entry and swap the other to the back.
        clusters(clustChoices(1)) = clusters(clustChoices(2)) + clusters(clustChoices(1));
        clusters([clustChoices(2) numNodes/2-q+1]) = clusters([numNodes/2-q+1 clustChoices(2)]);
        
        % Check if we've made a bigger cluster. If we have, update the
        % maxClust to reflect it.
        if clusters(clustChoices(1))/numNodes > maxClust(q)
            maxClust(q + 1) = clusters(clustChoices(1))/numNodes;
        else
            maxClust(q + 1) = maxClust(q);
        end
        
    end
    
end

% Calculate number of closed loops and maximum jump in the order parameter.
closedLoops(loopNum) = clusters(1);
maxJump = max(diff(maxClust));

% Now we do the part where we've only got closed loops left, which will eventually merge.
nodeProb = closedLoops/numNodes; % Probability of choosing a particular closed loop.
clProxy = closedLoops; % Proxy of closed loops for easier tracking.
while max(maxClust) < 1
    
    % Pull each node from a distribution scaled by size of the cluster.
    % Then remove each node chosen from the list.
    r = rand;
    nc1 = sum(r >= cumsum([0, nodeProb]));
    clProxy(nc1) = clProxy(nc1) - 1;
    nodeProb = clProxy/sum(clProxy);
    r = rand;
    nc2 = sum(r >= cumsum([0, nodeProb]));
    clProxy(nc2) = clProxy(nc2) - 1;
    nodeProb = clProxy/sum(clProxy);
    % Only do this if we're merging two closed loops.
    if nc1 ~= nc2
        closedLoops(nc1) = closedLoops(nc1) + closedLoops(nc2);
        maxClust(q + m) = max(maxClust(q + m - 1), closedLoops(nc1)/numNodes);
        closedLoops(nc2) = [];
        clProxy(nc1) = clProxy(nc1) + clProxy(nc2);
        clProxy(nc2) = [];
        m = m + 1;
        nodeProb = clProxy/sum(clProxy);
    else
        maxClust(q + m) = maxClust(q + m - 1);
        m = m + 1;
    end
    
    % Keep track of the maximum jump in the order parameter, as this is what we're interested in
    % finding.
    maxJump = max(diff(maxClust));
    
    % Figure out if we didn't finish connecting all the loops.
    if clProxy == 0
        nofinish = 1;
        break;
    end
end
end