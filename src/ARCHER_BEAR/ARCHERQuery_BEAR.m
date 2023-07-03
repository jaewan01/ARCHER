%
% ARCHERQuery_BEAR: Query phase of ARCHER using BEAR.
%
% Parameters
%   s : index of a seed node (use original index)
%   c : restarting probability
%   PREP : preprocessed matrices from BEAR
%   is_clique : if ARCHER selected clique exapnsion true, if ARCHER selected star expansion false
%   num_v : number of nodes
% Return values
%   r : RWR scores w.r.t. the seed node s (use original index)
%

function r = ARCHERQuery_BEAR(s, c, PREP, is_clique, num_v)

if is_clique
    %clique expansion based query phase
    r = BearQuery(s, c, PREP); 
else
    c_prime = 1 - sqrt(1-c);
    %star expansion based query phase
    rs = BearQuery(s, c_prime, PREP);
    r = c * rs(1:num_v, 1) / c_prime;
end

end

