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

function [PREP] = BePIPre(A, c, k)
% BePIPre: Pre-processing phase of BePI.
%
% Parameters
%   A : Adjacency matrix
%   c : restarting probability
%   k : hub selection ratio
%
% Return values
%   PREP: structure for preprocessing
%       invL1 : L1^{-1}
%       invU1 : U1^{-1}
%       S : schur complement
%       L2, U2: ilu(S)
%       H12, H21, H32, H32
%       order : reordered node indices (original index -> reordered index)
%

clear_flg = true;
nstep = 0;

[A, neworder, inv_idx, hubwidth, ndwidth] = SlashDeadBurn(A, nstep, k);

% row normalize adjacency matrix -> we give normalized matrix
st = tic;
%A = Normalize(A);

% compute H
n = size(A, 1);
H = (speye(n) - (1-c) * A');
%spy(H);

if clear_flg
    clearvars A;
end

n1 = ndwidth - hubwidth;
n2 = ndwidth;

H11 = H(1:n1, 1:n1);
H12 = H(1:n1, n1+1:n2);
H21 = H(n1+1:n2, 1:n1);
H22 = H(n1+1:n2, n1+1:n2);
H31 = H(n2+1:n, 1:n1);
H32 = H(n2+1:n, n1+1:n2);

if clear_flg
    clearvars H;
end
time_H = toc(st);

% LU decompose H_{11} and invert the results
st = tic;
[L1, U1] = lu(H11);
if clear_flg
    clearvars H11;
end

invL1 = inv(L1);
if clear_flg
    clearvars L1;
end

invU1 = inv(U1);
if clear_flg
    clearvars U1;
end
time_lu = toc(st);

% compute schur complement
st = tic;
S = H22 - H21*(invU1*(invL1*H12));
if clear_flg
    clearvars H22;
end 
time_S = toc(st);

setup.type = 'nofill';   

st = tic;
[L2, U2] = ilu(S, setup);
time_ilu = toc(st);

PREP.invL1 = invL1;
PREP.invU1 = invU1;
PREP.S = S;
PREP.L2 = L2;
PREP.U2 = U2; 
PREP.H12 = H12;
PREP.H21 = H21;
PREP.H31 = H31;
PREP.H32 = H32;
PREP.order = inv_idx;
end


function [Ak, new_idx, inv_idx, hubwidth, ndwidth] = SlashDeadBurn(A, nstep, ratio)
n = size(A, 1);
[Ak, new_idx, inv_idx, ndwidth] = Deadend(A, nstep);

[~, sborder, hubwidth] = SlashBurn(Ak(1:ndwidth, 1:ndwidth), ratio);

new_idx = new_idx([sborder',ndwidth+1:n]');

Ak = A(new_idx, new_idx);

inv_idx = zeros(n, 1);
inv_idx(new_idx) = (1:n)';
end

