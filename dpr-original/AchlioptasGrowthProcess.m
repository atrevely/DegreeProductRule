function [rEdge, cEdge, degs] = AchlioptasGrowthProcess(numNodes, numChoices, len, alpha, degDistPts)
% Performs an Achlioptas network growth process, for which edges compete
% for addition based on the product of cluster sizes. Used to recreate the
% results of Achlioptas et. al. in order to build a new growth process.

% numNodes = number of nodes in the network.
% numChoices = number of edges competing at each timestep.
% len = length of the run, in timesteps.
% alpha = parameter for interpolating between Erdos-Renyi and DPR.
% degDistPts = points for calculating the degree distribution.

% rEdge = first node of added edges, ordered by timestep.
% cEdge = second node of added edges, ordered by timestep.
% An order pair [rEdge(n),cEdge(n)] will specify an edge that has been
% added between two nodes at timestep n.
% degs = final degrees of each node.

% This version updates clusters manually rather than use scomponents, which
% makes is generally faster than Achlioptas version 1.

warning('off','all');
rng('shuffle');

% Initialize the sparse matrix of node connections, as well as degree
% vector and edge order vectors.
A = sparse(numNodes,numNodes);
choiceClstProd = zeros(1,numChoices);
rcNodeChoices = zeros(numChoices,2);
rEdge = zeros(1,len);
cEdge = zeros(1,len);
nodeCluster = scomponents(A);
degree = ones(1,numNodes);
degs = [];
m = 1;

% Network growth. Each iteration is a timestep that adds an edge.
for a = 1:len
    nC = 1;
    while nC <= numChoices
        % Samples in pairs with no replacement, meaning no self-loops.
        rcNodeChoices(nC,:) = sort(randperm(numNodes,2));
        % Check if edge is already occupied. If it is, then don't increment
        % counter so as to overwrite the existing edge choice.
        if A(rcNodeChoices(nC,1),rcNodeChoices(nC,2)) == 0 && length(unique(rcNodeChoices, 'rows')) == length(rcNodeChoices)
            nC = nC + 1;
        end
    end
    
    % Calculate competition criteria and choose the winning edge.
    for b = 1:numChoices
        choiceClstProd(b) = nnz(nodeCluster(rcNodeChoices(b,1))==nodeCluster).*nnz(nodeCluster(rcNodeChoices(b,2))==nodeCluster);
    end
    [~, locMin] = min(choiceClstProd);
    kappa = sum(choiceClstProd);
    probAdj = exp((min(choiceClstProd)-max(choiceClstProd))/(kappa*alpha));
    if rand(1) < probAdj
        locMin = randi(nC - 1);
    end
    
    % Add the winning edge to the edge tracking vectors and keep track of
    % the largest cluster currently in the network.
    rEdge(a) = rcNodeChoices(locMin,1);
    cEdge(a) = rcNodeChoices(locMin,2);
    A(rEdge(a), cEdge(a)) = 1;
    A(cEdge(a), rEdge(a)) = 1;
    if (nodeCluster(rEdge(a)) ~= nodeCluster(cEdge(a)))
        nodeCluster(nodeCluster == nodeCluster(cEdge(a))) = nodeCluster(rEdge(a));
    end
    
    % Bump up the degrees on the two nodes.
    degree(rcNodeChoices(locMin,1)) = degree(rcNodeChoices(locMin,1)) + 1;
    degree(rcNodeChoices(locMin,2)) = degree(rcNodeChoices(locMin,2)) + 1;
    
    % If we're at a point where we want the node degrees for creation of
    % the degree distribution, then record them here.
    if any(a == degDistPts)
        degs(m,:) = degree;
        m = m + 1;
    end
end
warning('on','all');
end