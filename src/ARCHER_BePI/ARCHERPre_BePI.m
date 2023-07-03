%
% ArcherPre_BePI: Pre-processing phase of ARCHER using BePI.
%
% Parameters
%   R : node-weight matrix
%   W : hyperedge-weight matrix
%   c : restart probability
% Return values
%   PREP : preprocessed matrices from BePI
%   is_clique : if ARCHER selected clique exapnsion true, if ARCHER selected star expansion false
%

function [PREP,is_clique] = ARCHERPre_BePI(R,W,c)
    % number of nodes
    num_v = size(W, 1);
    % number of hyperedges
    num_e = size(W, 2);

    %row_normalize W
    vecw = sum(W, 2);
    vecw = 1 ./ vecw;
    Dv_inv = spdiags(vecw, 0, num_v, num_v);
    W_norm = Dv_inv * W;

    %row_normalize R
    vecr = sum(R, 2);
    vecr = 1 ./ vecr;
    De_inv = spdiags(vecr, 0, num_e, num_e);
    R_norm = De_inv * R;
    
    %transition matrix for clique expansion
    P = W_norm * R_norm;
    
    %compute nonzero of H matrix for each method
    clique_nonzeros = nnz(P);
    star_nonzeros = nnz(W_norm) * 2 + num_v + num_e;
    
    %automatic selection part of ARCHER
    is_clique = false;
    if clique_nonzeros < star_nonzeros
        is_clique = true;
    end
    
    %hub selection ratio
    k = 0.2;
    
    if is_clique
        %clique expansion based preprocessing phase
        [PREP]=BePIPre(sparse(P), c, k);
    else
        %transition matrix for star expansion
        S = [sparse(num_v, num_v) W_norm];
        S = [S; [R_norm sparse(num_e, num_e)]];
        c_prime = 1 - sqrt(1-c);
        %star expansion based preprocessing phase
        [PREP]=BePIPre(sparse(S), c_prime, k);
    end

end
