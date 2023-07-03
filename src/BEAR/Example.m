%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    BEAR - Block Elimination Approach for Random Walk with Restart on Large Graphs.
%    Author: Anonymized
%    
%    Version: 1.0
%    Date: August 13, 2014
%
%    This software is free of charge under research purposes.
%    For commercial purposes, please contact the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This example shows how to run BEAR.
%%% It first loads an example graph, pre-process it, and compute RWR scores for 100 random seed nodes 

%%% Load an example graph and save it in the adjacency matrix format
data = csvread('Data.csv');
data = data + 1;
n = max(max(data));
m = size(data, 1);
A = sparse(data(:,1), data(:, 2), 1, n, n);

%%% Run the preprocess phase BEAR with restart probability c = 0.05 and drop tolerance = 0
c = 0.05;
tic;
[invL1, invU1, invL2, invU2, H12, H21, orderInfo]=BearPre(A, c, 0);
time = toc;

%%% print preprocessing time
fprintf('%f\n', time);

%%% computes RWR scores for 100 random query nodes 
n_q = 100;
query_nodes = randi(n, n_q, 1);
tic;
for q=1:n_q
    r = BearQuery(query_nodes(q), c, invL1, invU1, invL2, invU2, H12, H21, orderInfo); 
end
time = toc;

%%% print query time
fprintf('%f\n', time/n_q);


