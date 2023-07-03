%
% ARCHERQuery_BePI: Query phase of ARCHER using BePI.
%
% Parameters
%   s : index of a seed node (use original index)
%   c : restarting probability
%   epsilon: error tolerance
%   PREP : preprocessed matrices from BEAR
%   is_clique : if ARCHER selected clique exapnsion true, if ARCHER selected star expansion false
%   num_v : number of nodes
% Return values
%   r : RWR scores w.r.t. the seed node s (use original index)
%   iter: # of iterations
%

function [r, iter] = ARCHERQuery_BePI(s, c, epsilon, PREP, is_clique, num_v)

if is_clique
    %clique expansion based query phase
    [r, iter] = BePIQuery(s, c, epsilon, PREP); 
else
    c_prime = 1 - sqrt(1-c);
    %star expansion based query phase
    [rs, iter] = BePIQuery(s, c_prime, epsilon, PREP);
    r = c * rs(1:num_v, 1) / c_prime;
end

end
