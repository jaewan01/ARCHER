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
% BearPre: Pre-processing phase of BEAR.
%
% Parameters
%   A : Adjacency matrix
%   c : restarting probability
%   dropTol : drop tolerance
% Return values
%   invL1 : L_{1}^{-1}
%   invU1 : U_{1}^{-1}
%   invL2 : L_{2}^{-1} 
%   invU2 : U_{2}^{-1}
%   H12 : H_{12}
%   H21 : H_{21}
%   ordering : reordered node indices (original index -> reordered index)
%

function [PREP]=BearPre(A, c, dropTol) 
    % number of nodes
    n = size(A,1);

    % number of hubs removed at a time in SlashBurn
    k=floor(0.001*n);
    if k == 0
        k = 1;
    end

    % symmetrize A if it is an undirected network
    %{
    if sum(sum(triu(A, 1))) == 0 || sum(sum(tril(A, -1))) == 0
        A = 1.0 * ((A + A') > 0);
    end
    %}
    
    % reorder adjacency matrix using SlashBurn
    % SlashBurn was modified to reorder each diagonal block in the descending order of degrees 
    % A : reordered adjacency matrix
    % SBorder : reordered node indices by SlashBurn (reordered index -> original index)
    % wing_width : number of hubs
    [A,SBorder, wing_width] = SlashBurn(A, k);
    % rotate the result of slashburn 180 degrees so that hubs are located
    % after spokes and each diagonal block is reordered in the ascending
    % order of degrees
    A = rot90(A,2);
    
    % row normalize adjacency matrix
    %{
    vec = sum(A~=0, 2);
    vec = bsxfun(@max, vec, 1);
    vec = 1 ./ vec;
    n = length(vec);
    D = spdiags(vec(:),0,n,n);
    A = D * A;
    %}
    
    % compute H
    H = (speye(n) - (1-c) * A');
   
    % last index of spokes
    sparse_end = (n-wing_width);
    
    % first index of hubs
    dense_start = sparse_end+1;
    
    H11 = H(1:sparse_end, 1:sparse_end);
    H12 = H(1:sparse_end, dense_start:n);
    H21 = H(dense_start:n, 1:sparse_end);
    H22 = H(dense_start:n, dense_start:n);

    % LU decompose H_{11} and invert the results
    [L1, U1] = lu(H11);
    invL1 = inv(L1);
    invU1 = inv(U1);
    S = H22 - H21*(invU1*(invL1*H12));

    % Reorder S
    D = sum(S~=0,1)' + sum(S~=0,2);
    [~, Sorder] = sort(D);
    S = S(Sorder, Sorder);
    H12 = H12(:, Sorder);
    H21 = H21(Sorder, :);
    % LU decompose S and invert the results
    [L2, U2] = lu(S);
    invL2 = inv(L2);
    invU2 = inv(U2);
    
    % drop near-zero entries
    if dropTol > 0
        invL1 = invL1 .* (invL1 > dropTol) + invL1 .* (invL1 < - dropTol);
        invU1 = invU1 .* (invU1 > dropTol) + invU1 .* (invU1 < - dropTol);
        invL2 = invL2 .* (invL2 > dropTol) + invL2 .* (invL2 < - dropTol);
        invU2 = invU2 .* (invU2 > dropTol) + invU2 .* (invU2 < - dropTol);
        H12 = H12 .* (H12 > dropTol) + H12 .* (H12 < - dropTol);
        H21 = H21 .* (H21 > dropTol) + H21 .* (H21 < - dropTol);
    end
    
    % reordered index -> original index in H_{11}
    Horder = zeros(n,1);
    for i=1:n
        Horder((n+1)-i)=SBorder(i);
    end

    % combine reordered index in H_{11} and that in S
    Torder = Horder;
    base = dense_start - 1;
    for i=1:wing_width
         Torder(base+i) = Horder(base+Sorder(i));
    end
    
    % original index -> reordered index
    ordering = zeros(n,1);
    for i=1:n
        ordering(Torder(i))=i;
    end
    PREP.invL1 = invL1;
    PREP.invU1 = invU1;
    PREP.invL2 = invL2;
    PREP.invU2 = invU2;
    PREP.H12 = H12;
    PREP.H21 = H21;
    PREP.ordering = ordering;
end