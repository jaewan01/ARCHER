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
% SlashBurn: shatter graph and reorder nodes to make a compact adjacency matrix.
% Modified version for BEAR
%
% Parameter
%   AOrig : adjacency matrix of a graph.
%   k : # of nodes to to shatter in each iteration
%
% Return values
%   Ak : newly reordered adjacency matrix
%   newind : original index -> reordered index
%   wingwidth : number of hubs
%

function [Ak,newind, wingwidth] = SlashBurn(AOrig, k)

if nargin<2, error('SlashBurn requires at least 2 arguments.'); end

if ~ismatrix(AOrig)
	error('In SlashBurn, @B should be a matrix.');
end
if ~isscalar(k)
	error('In SlashBurn, @k should be a scalar integer.');
end
if k~=round(k)
	if k > 0 && k < 1
		k = round(size(AOrig,2)*k);
	else
		error('In SlashBurn, @k should be a scalar integer.');
	end
end

ASymm = AOrig | AOrig';
ASymm(ASymm>0) = 1;

niter=0;
n = size(AOrig,2);
totalind = zeros(1,n);
cur_lpos = 1;
cur_rpos = n;
gccind = 1:n;

gccsize = n;

while niter == 0 || gccsize > k
	niter = niter + 1;

	A = ASymm(gccind,gccind);
	[disind,newgccind,topind] = RemoveHighDegree(A, k);

	topind_size = length(topind);
	disind_size = length(disind);

	totalind(cur_lpos:cur_lpos + topind_size - 1) = gccind(topind);
	cur_lpos = cur_lpos + topind_size;
	totalind(cur_rpos - disind_size + 1:cur_rpos) = gccind(disind);
	cur_rpos = cur_rpos - disind_size;
	newind = totalind;

	gccind = gccind(newgccind);
	newind(cur_lpos:cur_rpos) = gccind;
	
	gccsize = size(gccind,2);
    
end

wingwidth = cur_lpos - 1;

Ak = AOrig(newind,newind);










