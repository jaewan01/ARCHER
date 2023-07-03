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


function [ A, n ] = Import( path, type, base, fid )
% Import : load the data graph
%
% Parameters
%   path : a path of the data
%   base : node index base of the graph (0 or 1)
%
% Return values
%   A : adjacency matrix of the graph
%   n : # of nodes
%

if nargin < 2
    type = 'csv';
end

if nargin < 3
    base = 0;
end

if nargin < 4
    fid = 1;
end

start_id = 1 - base;

if strcmp(type, 'csv')
    X = csvread(path);
elseif strcmp(type, 'tsv')
    X = dlmread(path);
end

X(:, 1) = X(:, 1) + start_id;
X(:, 2) = X(:, 2) + start_id;

if size(X, 2) < 3
   X = [X, ones(size(X, 1), 1)];
end

n = max(max( X(:, 1), X(:, 2) ));

A = sparse( X(:, 1), X(:, 2), X(:, 3), n, n);
A = A - spdiags(diag(A), 0, n, n); %remove self loops

end

