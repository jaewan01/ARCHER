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


function [ bs ] = GetSize( A )
%GETSIZE returns the size of the matrix A in MB
%
% Parameters
%   A : a matrix
%
% Return values
%   bs : the size of the matrix A in MB
%


s = whos('A');
bs = s.bytes/1024/1024;

end

