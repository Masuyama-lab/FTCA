% These are the codes of Fast Topological CIM-based Adaptive Resonance Theory (FTCA)
% proposed in "N. Masuyama, C. K. Loo, H. Ishibuchi, N. Amako, Y. Nojima, and Y. Liu, 
% "Fast topological adaptive resonance theory based on correntropy induced metric," 
% Proc. of 2019 IEEE Symposium Series on Computational Intelligence (SSCI 2019), 
% pp. 2215-2221, Xiamen, China, December 6-9, 2019."
% 
% Please contact "masuyama@omu.ac.jp" if you have any problems.
%  

% clc
clear
whitebg('black')

MIter = 1;     % Number of iterations

NR = 0.1; % Noise Rate [0-1]


load 2D_ClusteringDATASET; nData = 60000; 
DATA = [data(1:end,1) data(1:end,2)];
originDATA = DATA;

% Randamize data
ran = randperm(size(DATA,1));
DATA = DATA(ran,:);


% Normalization [0-1]
for k=1:size(DATA,2)
    mmin = min(DATA(:,k));
    mmax = max(DATA(:,k));
    DATA(:,k) = (DATA(:,k)-mmin) ./ (mmax-mmin);
end
DATA(isnan(DATA))=0;



% Parameters of FTCA =================================================
FTCAnet.numNodes    = 0;   % Number of clusters
FTCAnet.weight      = [];  % Mean of cluster
FTCAnet.CountNode = [];    % Counter for each node
FTCAnet.edge = zeros(2,2); % Initial connections (edges) matrix
FTCAnet.NewEdgedNode = []; % Node which creates new edge.
FTCAnet.adaptiveSig = [];

FTCAnet.minCIM = 0.2;
FTCAnet.Lambda = 50;      % Interval for Node deletion and topology construction
% ====================================================================




time_ftca=0;

for nitr = 1:MIter
    
    fprintf('Iterations: %d/%d\n',nitr,MIter);
    
    
    % Noise Setting [0,1]
    if NR > 0
        noiseDATA = rand(nData*NR, size(DATA,2));
        DATA(1:size(noiseDATA,1),:) = noiseDATA;
    end
    
    % Randamize data
    ran = randperm(size(DATA,1));
    DATA = DATA(ran,:);
    
    % FTCA ==============================================
    tic
    FTCAnet = FTCA(DATA, FTCAnet);
    time_ftca = time_ftca + toc;
    
    figure(1);
    myPlot(DATA, FTCAnet, 'FTCA');
    drawnow
    % ===================================================
    
    
end

