function [rEdge, cEdge, degs] = DPRGrowthProcess(numNodes, numChoices, len, alpha, degDistPts)
% Performs a network growth process realization for the Degree Product Rule Process,
% with a parameter alpha that allows for interpolating between Erdos-Renyi and the
% DPR Process. Outputs two ordered lists of the two nodes that are connected at
% each timestep, as well as the final degree of each node.

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

rng('shuffle');

% Initialize the degree vector, as well as the vectors holding the list of
% edges and the order they get added to the network.
degree = ones(1,numNodes);
rcNodeChoices = zeros(numChoices,2);
rEdge = zeros(1,len);
cEdge = zeros(1,len);
degs = [];
degDistCts = 1;

% Create a cell object of sparse arrays to keep track of edges.
edgeTrack = cell(numNodes,1);
atemp = cell(1);
atemp{1} = sparse(numNodes,1);
edgeTrack(:) = atemp;

% Network growth. Each iteration is a timestep that adds an edge.
for edgeNum = 1:len
    nC = 1; % Track how many of the requisite number of edge candidates we've selected so far.
    while (nC <= numChoices) && (nC <= (numNodes^2 - numNodes)/2 - edgeNum + 1)
        % Samples in pairs with no replacement, meaning no self-loops.
        rcNodeChoices(nC,:) = sort(randperm(numNodes,2));
        % Check if edge is already occupied. If it is, then don't increment
        % counter that we overwrite the existing edge choice until an unnocupied edge is selected.
        if nnz(edgeTrack{rcNodeChoices(nC,1)}(rcNodeChoices(nC,2))) == 0
            edgeTrack{rcNodeChoices(nC,1)}(rcNodeChoices(nC,2)) = 1;
            nC = nC + 1;
        end
    end
    choiceDegProd = degree(rcNodeChoices(:,1)).*degree(rcNodeChoices(:,2)); % Criteria for selecting the edge that gets added to the network.
    [~, locMin] = min(choiceDegProd);
    
    % Calculate the probability of switching our edge selection based on
    % the alpha input. alhpa = 0 no effect, alpha = inf Erdos-Renyi. Every
    % other value 0 < alpha < inf is something in between the two.
    kappa = sum(choiceDegProd);
    probAdj = exp((min(choiceDegProd)-max(choiceDegProd))/(kappa*alpha));
    if rand(1) < probAdj
        locMin = randi(nC - 1);
    end
    rEdge(edgeNum) = rcNodeChoices(locMin,1);
    cEdge(edgeNum) = rcNodeChoices(locMin,2);
    
    % Bump up the degrees on the two nodes.
    degree(rcNodeChoices(locMin,1)) = degree(rcNodeChoices(locMin,1)) + 1;
    degree(rcNodeChoices(locMin,2)) = degree(rcNodeChoices(locMin,2)) + 1;
    
    % Remove the edges not chosen from our edge tracker.
    rcNodeChoices(locMin,:) = [];
    for b = 1:nC-2
        edgeTrack{rcNodeChoices(b,1)}(rcNodeChoices(b,2)) = 0;
    end
    
    % If we're at a point where we want the node degrees for creation of
    % the degree distribution, then record them here.
    if any(edgeNum == degDistPts)
        degs(degDistCts,:) = degree;
        degDistCts = degDistCts + 1;
    end
end
end