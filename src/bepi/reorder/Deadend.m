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


%   This software is free of charge under research purposes.
%   For commercial purposes, please contact the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ Ak, new_idx, inv_idx, width ] = Deadend( AOrig, step )
% deadend performs deadend reordering
%
% Parameters
%   AOrig : Original adjacency matrix
%   step: step for subdangling node
%   
% Return values
%   Ak : Reordered matrix
%   new_idx: permutation for deadend reordering
%   inv_idx: inverse permuation of deadend reordering
%   width: # of deadends

if nargin < 2
   step = size(AOrig, 1); 
end

n = size(AOrig, 1);

nstep = 0;
nd_idx = (1:n)';

cur_rpos = n;
totalind = zeros(1,n);
new_idx = nd_idx;

while 1
    %spy(AOrig(new_idx, new_idx));
    nstep = nstep + 1;
    A = AOrig(nd_idx, nd_idx);
    D = sum(A~=0, 2);
    dead_idx = find(D == 0);
    nodead_idx = setdiff((1:length(nd_idx))', dead_idx);
    
    if isempty(dead_idx)
       break; 
    end
    
    dead_size = length(dead_idx);
    totalind(cur_rpos - dead_size + 1:cur_rpos) = nd_idx(dead_idx);
    cur_rpos = cur_rpos - dead_size;
    new_idx = totalind;
    
    nd_idx = nd_idx(nodead_idx);
    new_idx(1:cur_rpos) = nd_idx;
    
    if nstep > step
       break; 
    end
end

width = length(nd_idx);
Ak = AOrig(new_idx, new_idx);

inv_idx = zeros(n, 1);
inv_idx(new_idx) = (1:n)'; 


end

