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


function [Ak,newind, hubwidth, blksz] = SlashBurn(AOrig, ratio)
% SlashBurn: shatter graph and reorder nodes to make a compact adjacency matrix.
%
% Parameter
%   AOrig : adjacency matrix of a graph.
%   k : # of nodes to to shatter in each iteration
%
% Return values
%   Ak : newly reordered adjacency matrix
%   newind : original index -> reordered index
%   hubwidth : number of hubs
%   blksz : array of size of each block


if nargin<2, error('SlashBurn requires at least 2 arguments.'); end

if ~ismatrix(AOrig)
    error('In SlashBurn, @B should be a matrix.');
end
if ~isscalar(ratio)
    error('In SlashBurn, @k should be a scalar integer.');
end
if ratio~=round(ratio)
    if ratio > 0 && ratio < 1
        k = round(size(AOrig,2)*ratio);
    else
        error('In SlashBurn, @k should be a scalar integer.');
    end
end

ASymm = AOrig | AOrig';
%ASymm(ASymm>0) = 1;

niter=0;
n = size(AOrig,2);
totalind = zeros(1,n);
cur_lpos = 1;
cur_rpos = n;

gccind = 1:n;
gccsize = n;

%% added for block size
totblksize = zeros(1, n);
b_cur_rpos = n;

while niter == 0 || gccsize > k
   niter = niter + 1;
    
    A = ASymm(gccind,gccind);
    [disind,newgccind,topind, blksz] = RemoveHighDegree(A, k);
    
    topind_size = length(topind);
    disind_size = length(disind);    
    
    blksz_size = length(blksz); % added for block size
    totblksize(b_cur_rpos - blksz_size + 1:b_cur_rpos) = blksz; % added for block size
    b_cur_rpos = b_cur_rpos - blksz_size; % added for block size
    
    totalind(cur_lpos:cur_lpos + topind_size - 1) = gccind(topind);
    cur_lpos = cur_lpos + topind_size;
    totalind(cur_rpos - disind_size + 1:cur_rpos) = gccind(disind);
    cur_rpos = cur_rpos - disind_size;
    newind = totalind;
    
    gccind = gccind(newgccind);
    newind(cur_lpos:cur_rpos) = gccind;
    
    gccsize = size(gccind,2);
    k = ceil(ratio*gccsize);
    %if k < 10
    %    break;
    %end
end

hubwidth = cur_lpos - 1;
totblksize(b_cur_rpos) = gccsize;
blksz = totblksize(b_cur_rpos:end);

% to make hub section located in bottom-right area
newind = fliplr(newind);
blksz = fliplr(blksz);

if size(newind, 2) ~= 1
   newind = newind'; 
end

if size(blksz, 2) ~= 1
   blksz = blksz'; 
end

Ak = AOrig(newind,newind);










