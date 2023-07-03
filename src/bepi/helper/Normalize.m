%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BePI: Fast and Memory-Efficient Method for Billion-Scale Random Walk with Restart
%   
%   Authors: Jinhong Jung (jinhongjung@snu.ac.kr), Seoul National University
%            Namyong Park (namyong.park@snu.ac.kr), Seoul National University
%            Lee Sael (sael@sunykorea.ac.kr), The State University of New York (SUNY) Korea
%            U Kang (ukang@snu.ac.kr), Seoul National University
%   
%   Version : 1.1
%   Date: 2017-02-05
%   Main contact: Jinhong Jung
%
%   This code is released only for research purposes in universities.
%   To use this code for commericial purposes, contact U Kang (ukang@snu.ac.kr).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A, D]=Normalize(A) 
% Normalize: row-normalize given adjacency matrix
%
% Parameters
%   A : adjacency matrix
% Return values
%   A : normalized Adjacency matrix
%   D : diagonal matrix of degrees

    %vec = sum(A~=0, 2);
    vec = sum(A, 2);
    n = length(vec);  
    %D = spdiags(vec(:), 0, n, n);
    D = vec(:)';
    vec = bsxfun(@max, vec, 1);
    vec = 1 ./ vec;    
    invD = spdiags(vec(:),0,n,n);
    A = invD * A;

end
