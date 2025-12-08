function x = polar_encoder(info_bits, Q)
% info_bits: 1xK_total row (payload + CRC)
% Q: reliability ordering 1..N (length N)
N = length(Q);
K = length(info_bits);
n = log2(N);

info_pos = Q(N-K+1:end);
u = zeros(1,N);
u(info_pos) = info_bits;

% butterfly (F^{⊗n}) to get x_nat
x_nat = u;
m = 1;
for stage = 1:n
    for i = 1:2*m:N
        a = x_nat(i:i+m-1);
        b = x_nat(i+m:i+2*m-1);
        x_nat(i:i+2*m-1) = [mod(a+b,2), b];
    end
    m = 2*m;
end

% bit-reversal mapping
idx = zeros(1,N);
for k=0:N-1
    rev=0;
    for b=1:n
        rev = rev*2 + bitget(k,b);
    end
    idx(k+1)=rev+1;
end
x = x_nat(idx);
end
