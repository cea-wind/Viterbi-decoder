
trel = poly2trellis(7, [171 133]);
msg = randi([0,1],1,100);
code2 = convenc(msg,trel);

code = 1-2*code1 + 0.5*randn(size(code1));
tblen = 10;
decoded1 = vitbiDecoder(code,trel,tblen);

sum(abs(msg(1:end-tblen)-decoded1(tblen+1:end)))