function [ decoded ] = vitbiDecoder( code,trel,tblen )
%The simplest implementation of viterbi decoder
%This function can be optimazed in two aspect.
%   Using C program
%   Adjust algorithm by the structure of trellis
%The BER is lower of this function than vitdec in most cases. The reasons may includes
%   We choose the maximum metric at last.
%   Traceback depth largest.
K = log2(trel.numInputSymbols);
N = log2(trel.numOutputSymbols);
code = code(:).';
codeReshape = reshape(code,N,[])';
savedStates= zeros(trel.numStates,size(codeReshape,1)+1);
statesFlag = zeros(trel.numStates,size(codeReshape,1)+1);
decoded =  zeros(trel.numStates,size(codeReshape,1));
metric = zeros(trel.numStates,trel.numStates + 1);
% initial state is 0, some states can't reach if m < log2(trel.numStates)/K
statesFlag(1,1) = 1;
for m = 1:size(codeReshape,1)
	for n = 1:trel.numStates
		if(statesFlag(n,m)>0)
			for c = 1:trel.numInputSymbols
				selected = trel.nextStates(n,c);
				outputs = trel.outputs(n,c);
                outputs = dec2bin(outputs,N) - 48;
				metric_add = codeReshape(m,:)*(1-2*outputs)';
				if(statesFlag(selected+1,m+1)==0)
					metric(selected+1,m+1) = metric(n,m)+metric_add;
					decoded(selected+1,m) = c-1;
					statesFlag(selected+1,m+1) = n;
				else
					if(metric(n,m)+metric_add>metric(selected+1,m+1))
						metric(selected+1,m+1) = metric(n,m)+metric_add;
                        decoded(selected+1,m) = c-1;
                        statesFlag(selected+1,m+1) = n;
					end
				end
			end
		end
	end
end

[~,m] = max(metric(:,end));
decoded_seq = [];
for n=size(codeReshape,1):-1:1
    decoded_seq = [dec2bin(decoded(m,n),K)-48 decoded_seq];
    m = statesFlag(m,n+1);
end
decoded = [zeros(1,tblen) decoded_seq(1:end-tblen)];
end