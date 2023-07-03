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
% RemoveHighDegree: remove k hubs from GCC
% Modified version for BEAR
%
% Parameter
%   B : giant connected component (GCC) 
%   k : # of hubs removed at a time
%
% Return values
%   disind : indices of nodes in disconnected component
%   gccind : indices of nodes in GCC
%   topind : indices of top k nodes
%
% Example:
%
function [disind,gccind,topind] = RemoveHighDegree(B, k)
	
if k == 1
	[~,topk] = max(sum(B,2));
else 
	[~,I] = sort(sum(B,2), 'descend');
	topk = I(1:k);
end

topk = topk(:);

% Variable: Type and Size; Range of Value; Meaning
%	S: integer; > 0; # connected components (CCs)
%	C: integer vector of size n; 1~S; CC label assigned to each node
%	H: integer vector of size S; 0 >; # nodes having each label
%	T: integer vector of size S; 0~k; topk node id to which each CC is attached
[S, C, H, T] = ComputeConnComp(B, topk);

% Random suffling of component labels
% This needs because ComputeConnComp produces labeling undesirably ordered 
%	by top-k even when @dis_select == 'dis_naive'
randC = randperm(S);
H = H(randC);
T = T(randC);
randmap = zeros(S,1);
randmap(randC) = 1:S;
C = randmap(C);

[~,gcc] = max(H);
topind = topk;
gccind = find(C==gcc);
disind = find(C~=gcc);
disind = disind(~ismember(disind, topind));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified part of SlashBurn for BEAR
% Reorder each diagonal block in the descending order of degrees
% (It will be rotated later so that each diagonal block is reordered in the
% ascending order of degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degree_to_topK = sum(B(topind,:))' + sum(B(:, topind),2);

if ~isempty(disind)
	degree = sum(B,2);
	degree = degree + sum(B,1)';
    degree = degree - degree_to_topK;
	degree = degree(disind); %within block degree
    disind_info = [T(C(disind)) H(C(disind)) C(disind) degree];
	colspec = [-1 -2 3 -4];
	[~,disindI] = sortrows(disind_info, colspec);
	disind = disind(disindI);
end

end