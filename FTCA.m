% 
% (c) 2021 Naoki Masuyama
% 
% These are the codes of Fast Topological CIM-based Adaptive Resonance Theory (FTCA)
% proposed in "N. Masuyama, C. K. Loo, H. Ishibuchi, N. Amako, Y. Nojima, and Y. Liu, 
% "Fast topological adaptive resonance theory based on correntropy induced metric," 
% Proc. of 2019 IEEE Symposium Series on Computational Intelligence (SSCI 2019), 
% pp. 2215-2221, Xiamen, China, December 6-9, 2019."
% 
% Please contact "masuyama@omu.ac.jp" if you have any problems.
%  
function net = FTCA(DATA, net)

numNodes = net.numNodes;       % Number of nodes
weight = net.weight;           % Node position
CountNode = net.CountNode;     % Counter for each node
adaptiveSig = net.adaptiveSig; % sigma in each node

% Parameters for Topology
edge = net.edge;
NewEdgedNode = net.NewEdgedNode;
Lambda = net.Lambda;

minCIM = net.minCIM;


for sampleNum = 1:size(DATA,1)
    
    % Calculate an initial adaptiveSig based on inputs.
    if isempty(weight) == 1 && ((sampleNum-1) + Lambda <= size(DATA,1))
        % Extract inputs
        exNodes = DATA( sampleNum:(sampleNum-1) + Lambda , :);
        
        % Normalization [0-1]
%         for k=1:size(exNodes,2)
%             mmin = min(exNodes(:,k));
%             mmax = max(exNodes(:,k));
%             exNodes(:,k) = (exNodes(:,k)-mmin) ./ (mmax-mmin);
%         end
%         exNodes(isnan(exNodes))=0;
        
%         exNodes = zscore(exNodes);
        
        [n,d] = size(exNodes);
        exNodesStd = std(exNodes);
        
        % normal reference rule-of-thumb
        % https://www.sciencedirect.com/science/article/abs/pii/S0167715212002921
        estSig = median( ((4/(2+d))^(1/(4+d))) * exNodesStd * n^(-1/(4+d)) );
        
    elseif mod(sampleNum, Lambda) == 1 && (sampleNum - Lambda) >= 0
        % Extract inputs
        exNodes = DATA( (sampleNum+1) - Lambda:sampleNum, :);
        
        % Normalization [0-1]  
%         for k=1:size(exNodes,2)
%             mmin = min(exNodes(:,k));
%             mmax = max(exNodes(:,k));
%             exNodes(:,k) = (exNodes(:,k)-mmin) ./ (mmax-mmin);
%         end
%         exNodes(isnan(exNodes))=0;
        
%         exNodes = zscore(exNodes);
        
        [n,d] = size(exNodes);
        exNodesStd = std(exNodes);
        
        % normal reference rule-of-thumb
        % https://www.sciencedirect.com/science/article/abs/pii/S0167715212002921
        estSig = median( ((4/(2+d))^(1/(4+d))) * exNodesStd * n^(-1/(4+d)) );
        
    end
    
    
    
    
    % Current data sample.
    input = DATA(sampleNum,:);
    
    if size(weight,1) < 2 % In the case of the number of nodes in the entire space is small.
        % Add Node
        numNodes = numNodes + 1;
        weight(numNodes,:) = input;
        CountNode(numNodes) = 1;
        NewEdgedNode(1, numNodes) = 0;
        edge(numNodes, :) = 0;
        edge(:, numNodes) = 0;
        adaptiveSig(numNodes) = estSig;
        
    else
        
        % Calculate CIM based on global mean adaptiveSig.
        globalCIM = CIM(input, weight, mean(adaptiveSig));
        
        % Extract local area nodes which placed on around input.
        [~, idxLocal] = find(globalCIM < 1.0-minCIM);
        
        
        if size(idxLocal,2) < 2 % In the case of the number of nodes around input is small.
            % Add Node
            numNodes = numNodes + 1;
            weight(numNodes,:) = input;
            CountNode(numNodes) = 1;
            NewEdgedNode(1, numNodes) = 0;
            edge(numNodes,:) = 0;
            edge(:,numNodes) = 0;
            adaptiveSig(numNodes) = mean(adaptiveSig);
            
        else
            
            % Sorting by CIM state for finding the winner node.
            localCIM = globalCIM(idxLocal);
            [~, orderLocal] = sort(localCIM, 'ascend');
            
            % Set CIM state of 1st local area winner node.
            Lcim1 = localCIM(1, orderLocal(1));
            
            % For resonance process.
            resonance = false;
            currentSortedIndex = 1;
            
            
            if minCIM < Lcim1 % In the case of input is far from winner nodes.
                % Add Node
                numNodes = numNodes + 1;
                weight(numNodes,:) = input;
                CountNode(numNodes) = 1;
                NewEdgedNode(1, numNodes) = 0;
                edge(numNodes,:) = 0;
                edge(:,numNodes) = 0;
                adaptiveSig(numNodes) = min(adaptiveSig(idxLocal)); 
                
            else
                
                while ~resonance
                    
                    % Set indexes of local winner nodes s1 and s2 by currentSortedIndex.
                    s1 = idxLocal( orderLocal(currentSortedIndex) );
                    s2 = idxLocal( orderLocal(currentSortedIndex + 1) );
                    
                    % Set CIM state between the local winner nodes and the input for Vigilance Test.
                    Lcim_s1 = localCIM( 1, orderLocal(currentSortedIndex) );
                    Lcim_s2 = localCIM( 1, orderLocal(currentSortedIndex + 1) );
                    
                    % Vigilance Test
                    if Lcim_s1 < minCIM && Lcim_s2 < 1.0-minCIM
%                     if Lcim_s1 < minCIM
                        % Match Success
                        weight(s1,:) = weight(s1,:) + (1/CountNode(s1)) * (input - weight(s1,:));
                        CountNode(s1) = CountNode(s1) + 1;
                        
                        if Lcim_s2 < minCIM
                            % Update weight of s1 neighbors.
                            s1Neighbors = find(edge(s1,:));
                            for k = s1Neighbors
                                weight(k,:) = weight(k,:) + (1/(10*CountNode(k))) * (weight(s1,:) - weight(k,:));
                            end
                            
                            % Create an edge between s1 and s2 nodes.
                            NewEdgedNode(1,s1) = 1;
                            edge(s1,s2) = 1;
                            edge(s2,s1) = 1;
                            
                            % Sigma adaptation. (For a Journal Improvement.) ===============================
%                             if nnz(edge(s1,:)) ~= 0
%                                 n1 = CIM(weight(s1,:), weight(s2,:), adaptiveSig(s1));
%                                 n2 = CIM(weight(s2,:), weight(s1,:), adaptiveSig(s2));
%                                 if minCIM < n1 && minCIM < n2 % Consider between s1 and s2 nodes, not an input and nodes.
%                                     adaptiveSig(s1) = adaptiveSig(s1) - (Lcim_s1 * Lcim_s2) * (1/CountNode(s1)) * adaptiveSig(s1);
%                                 end
%                             end
                            % ==============================================================================
                            
                        end
                        
                        resonance = true;
                        
                    else
                        % Match Fail
                        if(currentSortedIndex == size(idxLocal,2)-1) % If searched all the existing nodes, then...
                            % Add Node
                            numNodes = numNodes + 1;
                            weight(numNodes,:) = input;
                            CountNode(numNodes) = 1;
                            NewEdgedNode(1, numNodes) = 0;
                            edge(numNodes, :) = 0;
                            edge(:, numNodes) = 0;
                            adaptiveSig(numNodes) = min(adaptiveSig(idxLocal));
                            
                            resonance = true;
                            
                        else
                            currentSortedIndex = currentSortedIndex + 1;    % Search another cluster orderd by sortedNodes.
                        end
                        
                    end % if c1 < minCIM && c2 < 1.0
                    
                end % while ~resonance
                
                
                % Plotting per one sample.
%                 if mod(sampleNum, 1) == 0
%                     % Delete Intersections of edge
%                     [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));
%                     
%                     net.numNodes = numNodes;      % Number of nodes
%                     net.weight = weight;          % Mean of nodes
%                     net.CountNode = CountNode;    % Counter for each node
%                     net.edge = edge;
%                     net.Lambda = Lambda;
%                     figure(3);
%                     myplot(DATA, net, 'TCA Euclidean');
%                     hold on
% %                     plot(input(:,1),input(:,2),'w+','Markersize',15, 'Linewidth',15);
%                     hold off
%                 end
                
                
            end % if minCIM < Lcim1
            
        end % if size(idxEx,2) < 2
        
    end % if size(weight,1) < 2
        
    
    
    % Topology Reconstruction
    if mod(sampleNum, Lambda) == 0
        % -----------------------------------------------------------------
        % Delete Node based on number of neighbors
        nNeighbor = sum(edge);
        deleteNodeEdge = (nNeighbor == 0);
        
        % Delete process
        numNodes = numNodes - sum(deleteNodeEdge);
        weight(deleteNodeEdge, :) = [];
        CountNode(deleteNodeEdge) = [];
        NewEdgedNode(:, deleteNodeEdge) = [];
        edge(deleteNodeEdge, :) = [];
        edge(:, deleteNodeEdge) = [];
        adaptiveSig(deleteNodeEdge) = [];
        
        % -----------------------------------------------------------------
        % Delete Intersections of edge
        [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));
        
    end % if mod(sampleNum, Lambda) == 0
    
    
end % for sampleNum = 1:size(DATA,1)



% -----------------------------------------------------------------
% Delete Node based on number of neighbors
nNeighbor = sum(edge);
deleteNodeEdge = (nNeighbor == 0);

% Delete process
numNodes = numNodes - sum(deleteNodeEdge);
weight(deleteNodeEdge, :) = [];
CountNode(deleteNodeEdge) = [];
NewEdgedNode(:, deleteNodeEdge) = [];
edge(deleteNodeEdge, :) = [];
edge(:, deleteNodeEdge) = [];
adaptiveSig(deleteNodeEdge) = [];

% -----------------------------------------------------------------
% Delete Intersections of edge
[weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));


% Cluster Labeling based on edge (Functions are available above R2015b.)
connection = graph(edge ~= 0);
LabelCluster = conncomp(connection);



net.numNodes = numNodes;      % Number of nodes
net.weight = weight;          % Mean of nodes
net.CountNode = CountNode;    % Counter for each node
net.adaptiveSig = adaptiveSig;
net.LabelCluster = LabelCluster;
net.edge = edge;
net.NewEdgedNode = NewEdgedNode;
net.Lambda = Lambda;

end


% Correntropy induced Metric (Gaussian Kernel based)
function cim = CIM(X,Y,sig)
% X : 1 x n
% Y : m x n
[n, att] = size(Y);
g_Kernel = zeros(n, att);

for i = 1:att
    g_Kernel(:,i) = GaussKernel(X(i)-Y(:,i), sig);
end

ret0 = GaussKernel(0, sig);
ret1 = mean(g_Kernel, 2);

cim = sqrt(ret0 - ret1)';
end

function g_kernel = GaussKernel(sub, sig)
% g_kernel = exp(-sub.^2/(2*sig^2));
g_kernel = 1/(sqrt(2*pi)*sig) * exp(-sub.^2/(2*sig^2));
end



% Delete intersections of edge
function [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, sigma)

% for d = 1:size(weight,1); % Search all nodes
for d = find(NewEdgedNode == 1) % Search only new edged nodes
    
    node1 = find(edge(d,:)); % Neighbors of d-th node
    if size(node1,1) >= 1
       posX1 = weight(d,:); % position of d-th node
        for m = 1:size(node1,2) % Search all neighbors of d-th nodes
            posY1 = weight(node1(m),:); % position of m-th neighbor node of d-th node
            for h = 1:size(node1,2)
                target2 = node1(h);
                node2 = find(edge(target2,:)); % Neighbors of m-th node
                posX2 = weight(target2,:); % position of h-th neighbor node of m-th node
                for k = 1:size(node2,2)
                    posY2 = weight(node2(k),:); % position of k-th neighbor node of h-th node
                    isConvex = findIntersection(posX1, posY1, posX2, posY2); % find intersections
                    if isConvex == 1 % If intersection is exist, delete edge which has larger CIM.
                        cim1 = CIM(weight(d,:), weight(node1(m),:), sigma);
                        cim2 = CIM(weight(target2,:), weight(node2(k),:), sigma);
                        if cim2 >= cim1
                            edge(target2, node2(k)) = 0;
                            edge(node2(k), target2) = 0;
                        else
                            edge(d, node1(m)) = 0;
                            edge(node1(m), d) = 0;
                        end
                    end % end isConvex
                end % end k
            end % end h
        end % end m  
    end

end % end d

NewEdgedNode = zeros(size(NewEdgedNode));

end

% Check intersection of edges
function [isConvex] = findIntersection(A, B, C, D)

F1  = B(:,1)-D(:,1);
F2  = B(:,2)-D(:,2);
M11 = B(:,1)-A(:,1);
M21 = B(:,2)-A(:,2);
M12 = C(:,1)-D(:,1);
M22 = C(:,2)-D(:,2);
deter = M11.*M22 - M12.*M21;
lambda = -(F2.*M12-F1.*M22)./deter;
gamma = (F2.*M11-F1.*M21)./deter;

% E = (lambda*[1 1]).*A + ((1-lambda)*[1 1]).*B;
% isConvex = (0 <= lambda & lambda <= 1)  & (0 <= gamma & gamma <= 1);

isConvex = (0 < lambda & lambda < 1)  & (0 < gamma & gamma < 1) ;
isConvex = isConvex';

end




