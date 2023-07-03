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

% 
% BearQuery: Query phase of BEAR.
%
% Parameters
%   s : index of a seed node (use original index)
%   c : restarting probability
%   invL1 : L_{1}^{-1}
%   invU1 : U_{1}^{-1}
%   invL2 : L_{2}^{-1} 
%   invU2 : U_{2}^{-1}
%   H12 : H_{12}
%   H21 : H_{21}
%   ordering : reordered node indices (original index -> reordered index)
%
% Return values
%   r : RWR scores w.r.t. the seed node s (use original index)
%

function [r]=BearQuery(s, c, PREP)
    
    invL1 = PREP.invL1;
    invU1 = PREP.invU1;
    invL2 = PREP.invL2;
    invU2 = PREP.invU2;
    H12 = PREP.H12;
    H21 = PREP.H21;
    ordering = PREP.ordering;

    % number of nodes
    n = size(invL1,1)+size(invL2,1);

    % last index of spokes
    sparse_end = size(invL1);

    % first index of hubs
    dense_start = sparse_end+1;

    % initialize starting vector
    q = sparse(n, 1);
    q(ordering(s))=1;

    % partition starting vector
    q1 = q(1:sparse_end, :);
    q2 = q(dense_start:n, :);

    % compute RWR score vector
    q_tilda = q2 - H21*(invU1*(invL1*q1));
    r2 = (invU2*(invL2*q_tilda));
    r1 = invU1*(invL1 *(q1 - H12*r2));
    r = c*[r1; r2];

    %reorder RWR scroe vector according to the original ordering
    r = full(r);
    r = r(ordering);

end


