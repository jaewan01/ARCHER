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


function [disind,gccind,topind,blksz] = RemoveHighDegree(B, k)
% RemoveHighDegree: remove k hubs from GCC
%
% Parameter
%   B : giant connected component (GCC) 
%   k : # of hubs removed at a time
%
% Return values
%   disind : indices of nodes in disconnected component
%   gccind : indices of nodes in GCC
%   topind : indices of top k nodes
%   blksz : array of block size
	
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
%randC = randperm(S);
%H = H(randC);
%T = T(randC);
%randmap = zeros(S,1);
%randmap(randC) = 1:S;
%C = randmap(C);

[~,gcc] = max(H);
topind = topk;
gccind = find(C==gcc);
disind = find(C~=gcc);
disind = disind(~ismember(disind, topind));

degree_to_topK = sum(B(topind,:))' + sum(B(:, topind),2);

blksz = [];

if ~isempty(disind)
	degree = sum(B,2);
	degree = degree + sum(B,1)';
    degree = degree - degree_to_topK;
	degree = degree(disind); %within block degree
    disind_info = [T(C(disind)) H(C(disind)) C(disind) degree];
	colspec = [-1 -2 3 -4];
	[~,disindI] = sortrows(disind_info, colspec);
	disind = disind(disindI);
    %blksz = H(unique(C(disind), 'stable')); % added for blksz << no supported in R2011b
	blksz = H(unique(C(disind))); % << should be modified in this part
end

end
