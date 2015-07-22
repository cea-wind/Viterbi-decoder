function code = convEncoder(msg,trel)
%The simplest encoding function of convolution codes.
%   WARNING : The running speed of convEncoder is slower than convenc.
%   because we didn't use c program to optimize it.
%   msg: vector, size 1*n
%   trel: the return structure of cctrellis or poly2trellis
%   code: vector, size 1*m
K = log2(trel.numInputSymbols);
N = log2(trel.numOutputSymbols);
msgReshape = reshape(msg,K,[])';
msgDec = bin2dec(num2str(msgReshape));
currentState = 0;
code = zeros(length(msgDec),N);
for n = 1:length(msgDec)
    code(n,:) = dec2bin(trel.outputs(currentState+1,msgDec(n)+1),N)-48;
    currentState = trel.nextStates(currentState+1,msgDec(n)+1);
end
code = reshape(code',1,[]);
end

