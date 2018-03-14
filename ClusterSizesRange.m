function [fracConn] = ClusterSizesRange(rEdge, cEdge, numNodes, sPoint, ePoint)
% Calculate the largest cluster size of a sliding window interval, with
% step size always one link. Uses manual method of finding clusters rather
% than scomponents at every step.

% numNodes = number of nodes in the network.
% numChoices = number of edges competing at each timestep.
% sPoint = timestep where we'd like to start calculating the order
% parameter.
% ePoint = timestep where we're done calculating the order parameter (and
% run stops).
% cType = type of network growth.
% alpha = paramete for interpolating between Erdos-Renyi and DPR.
% distPts = points for calculating the cluster size distribution.

% fracConn = fraction of nodes belonging to the largest cluster at each
% timestep.

% This works equivalently for forward and backward percolation.

warning('off','all');

% We want to start at sPoint by computing the clusters, then go
% manually from there forward to the end of the specified range. This is
% generally faster than recalculating the largest cluster by using
% scomponents from the gaimc package.
fracConn = zeros(1,ePoint - sPoint);
rEdgeThru = rEdge(1:sPoint);
cEdgeThru = cEdge(1:sPoint);
A = sparse(rEdgeThru,cEdgeThru,ones(1,sPoint),numNodes,numNodes);
A = spones(A);
A = round(0.5*(A + A'));
A(find(speye(numNodes))) = 0;
nodeCluster = scomponents(A);
largestClst = sum(nodeCluster == mode(nodeCluster))/numNodes;
fracConn(1) = largestClst;

% Now manually change the cluster assignments after each link is added.
for q = 2:(ePoint - sPoint)
    rnodeAdded = rEdge(sPoint + q - 1);
    cnodeAdded = cEdge(sPoint + q - 1);
    if (nodeCluster(rnodeAdded) ~= nodeCluster(cnodeAdded))
        nodeCluster(nodeCluster == nodeCluster(cnodeAdded)) = nodeCluster(rnodeAdded);
        largestClst = max(largestClst, length(nodeCluster(nodeCluster == nodeCluster(cnodeAdded))));
    end
    fracConn(q) = largestClst/numNodes;
end
warning('on','all');
end