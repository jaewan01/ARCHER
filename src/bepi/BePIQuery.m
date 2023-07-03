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


function [r, iter] = BePIQuery(s, c, epsilon, PREP)
% BePIQuery: Query phase of BePI.
%
% Parameters
%   s: seed node
%   c: restart probability
%   epsilon: error tolerance
%   PREP: structure for preprocessing
%       invL1 : L1^{-1}
%       invU1 : U1^{-1}
%       S : schur complement
%       L2, U2: ilu(S)
%       H12, H21, H32, H32
%       order : reordered node indices (original index -> reordered index)
%
% Return values
%   r: RWR score vector w.r.t. s
%   iter: # of iterations
%


nrestart = 5;

invL1 = PREP.invL1;
invU1 = PREP.invU1;
S = PREP.S;
L2 = PREP.L2;
U2 = PREP.U2;
H12 = PREP.H12;
H21 = PREP.H21;
H31 = PREP.H31;
H32 = PREP.H32;
order = PREP.order;

n1 = size(invL1, 1);
n2 = size(H21,1);
n3 = size(H31,1);
n = n1 + n2 + n3;

% initialize starting vector
L = length(s);

q = sparse(n, 1);
q(order(s))=1/L;

q1 = q(1:n1, :);
q2 = q(n1+1:n1+n2, :);
q3 = q(n1+n2+1:n, :);

q_tilda = (q2 - H21*(invU1*(invL1*q1)));
r2 = q_tilda;
iter = 0;
if nnz(q_tilda) ~= 0
    [r2, flag, relres, iter] = gmres(S, q_tilda, nrestart, epsilon, 100, L2, U2, q_tilda);
    %display(relres)
end
r1 = invU1*(invL1 *(q1 - H12*r2));
r3 = q3 - H31*r1 - H32*r2;
r = c*[r1;r2;r3];
r = full(r);
r = r(order);
r = r./sum(r); % for dead-end

if length(iter) > 1
    %iter = iter(1) * iter(2);
    iter = (iter(1)-1)*nrestart + iter(2);
    %iter = iter(2);
end

end
