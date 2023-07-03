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
% Mex_compile: Pre-processing phase of BEAR.
%

function Mex_compile(debug_mode)

if nargin == 0
	debug_mode = 0;
end

if debug_mode == 1
	debug_opt = ', ''-g''';
else
	debug_opt = '';
end


compile_str = {
	strcat('mex(''cpp/ComputeConnComp.cpp'', ''-largeArrayDims''', debug_opt,');'), ...
	strcat('mex(''cpp/ComputeTopkGreedy.cpp'', ''-largeArrayDims''', debug_opt,');')
};

for i=1:length(compile_str)
	fprintf('Evaluate:\n\t%s\n', compile_str{i});
	eval(compile_str{i});
	fprintf('Done!\n');
end