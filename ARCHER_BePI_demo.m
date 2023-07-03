W_raw = importdata('Dataset/sample_w.txt');
W = sparse(W_raw(:,1),W_raw(:,2),W_raw(:,3));

R_raw = importdata('Dataset/sample_r.txt');
R = sparse(R_raw(:,1),R_raw(:,2),R_raw(:,3));

edge_sizes = sum(W~=0,1);
mean_edge = mean(edge_sizes);
max_edge = max(edge_sizes);
fprintf("average edge size : %f, maximum edge size : %d\n",full(mean_edge), full(max_edge));

clearvars -except W R;

addpath('src/ARCHER_BePI');
addpath(genpath('src/bepi'));

c = 0.05;
k = 0.2;

%number of preprocessing times
n_p = 10;
%number of query times
n_q = 30;

prep_times = zeros(1, n_p);
query_times = zeros(1, n_q);

num_v = size(W, 1);
query_nodes = randi(num_v, n_q, 1);

fprintf("ARCHER_BePI Prep times: ");

for p=1:n_p
    clearvars -except W R c k n_p n_q prep_times query_times p query_nodes
    tic;
    [PREP, is_clique] = ARCHERPre_BePI(R,W,c);
    prep_time = toc;
    prep_times(1, p) = prep_time;
    fprintf("%.4f\t",prep_time);
end
fprintf("\n");

epsilon = 1e-9;

fprintf("ARCHER_BePI Query times: ");

num_v = size(W, 1);
for q=1:n_q
    tic;
    r = ARCHERQuery_BePI(query_nodes(q), c, epsilon, PREP, is_clique, num_v); 
    query_time = toc;
    query_times(1, q) = query_time;
    fprintf("%.4f\t",query_time);
end
fprintf("\n");

rmpath('src/ARCHER_BePI');
rmpath(genpath('src/bepi'));

clearvars -except prep_times query_times;

fprintf("ARCHER_BePI Prep mean:%.4f, std:%.4f\n", mean(prep_times(1,1:end)), std(prep_times(1,1:end)));
fprintf("ARCHER_BePI Query mean:%.4f, std:%.4f\n", mean(query_times(1,1:end)), std(query_times(1,1:end)));
