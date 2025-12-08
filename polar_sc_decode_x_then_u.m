function [u_hat, msg_hat] = polar_sc_decode_x_then_u(Lx, Q, K_total)
% Successive Cancellation (SC) decoder for Polar Codes
%
% Inputs:
%   Lx      : LLRs for encoded x in BIT-REVERSED order (1 x N)
%   Q       : Reliability sequence of length N (1..N)
%   K_total : Total info bits (payload + CRC)
%
% Outputs:
%   u_hat   : Decoded u-vector (1 x N) in NATURAL order
%   msg_hat : Extracted info bits in correct order
%
% NOTE:
%   This SC decoder follows the standard Arıkan f/g recursion
%   with iterative tree traversal.

N = length(Lx);
n = log2(N);

%% ---------- Undo Bit Reversal ----------
idx = zeros(1, N);
for k = 0:N-1
    rev = 0;
    for b = 1:n
        rev = rev*2 + bitget(k, b);
    end
    idx(k+1) = rev + 1;
end
L0 = Lx(idx);   % natural-order LLRs

%% ---------- Reliability: info + frozen ----------
info_pos   = Q(N-K_total+1 : end);
frozen_pos = Q(1 : N-K_total);

frozen_mask = false(1, N);
frozen_mask(frozen_pos) = true;

%% ---------- Allocate LLR and partial sum arrays ----------
LLR = cell(n+1, 1);
PS  = cell(n+1, 1);

for lev = 1:n+1
    LLR{lev} = zeros(1, N);
    PS{lev}  = zeros(1, N);
end
LLR{1} = L0;  % root node

u_hat = zeros(1, N);  % final decoded bits

%% ---------- Helper f and g ----------
f = @(a,b) sign(a).*sign(b).*min(abs(a),abs(b));
g = @(a,b,uhat) b + ((-1).^uhat).*a;

%% ---------- Main SC loop over each bit ----------
for i = 1:N

    %% --- Step 1: compute LLR for current bit by descending the tree ---
    level = 1;
    L_tmp = LLR{level};

    for lev = 2:n+1
        block = 2^(n-(lev-1));
        k = floor((i-1)/block);
        offset = k*block + 1;

        half = block/2;

        if mod(i-1, block) < half
            L_tmp = f( L_tmp(offset:offset+half-1), ...
                       L_tmp(offset+half:offset+block-1) );
        else
            u_left = PS{lev}(offset:offset+half-1);
            L_tmp = g( L_tmp(offset:offset+half-1), ...
                       L_tmp(offset+half:offset+block-1), ...
                       u_left );
        end
    end

    L_i = L_tmp(1);

    %% --- Step 2: hard decision ---
    if frozen_mask(i)
        u_hat(i) = 0;
    else
        u_hat(i) = (L_i < 0);
    end

    %% --- Step 3: update partial sums going upward ---
    PS{n+1}(i) = u_hat(i);

    for lev = n:-1:1
        block = 2^(n-(lev-1));
        k = floor((i-1)/block);
        offset = k*block + 1;
        half = block/2;

        left_range  = offset : offset+half-1;
        right_range = offset+half : offset+block-1;

        if mod(i-1, block) < half
            % update only right child PS via XOR
            PS{lev}(right_range) = PS{lev}(right_range) ...
                                 .^ PS{lev+1}(left_range);
        else
            PS{lev}(right_range) = PS{lev+1}(right_range);
        end
    end

end

%% ---------- Extract info bits ----------
msg_hat = u_hat(info_pos);

end
