G = [1 1 1;1 0 1];
[N,L] = size(G);
M = L-1;Ns = 2^M;
for state_i = 1:Ns
    state_b = dec2bin(state_i-1,M)-48;
    for input_bit = 0:1
        d_k = input_bit;
        a_k = rem(G(1,:)*[d_k state_b]',2);
        out(input_bit+1,:) = [d_k rem(G(2,:)*[a_k state_b]',2)];
        state(input_bit+1,:) = [ak state_b(1:M-1)];
    end
    nout(state_i,:) = 2*[out(1,:) out(2,:)]-1;
end