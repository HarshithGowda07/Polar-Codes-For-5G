function x = polar_encode_butterfly(info_bits)
% POLAR_ENCODE_BUTTERFLY  Polar encoder using butterfly structure (no matrix).
% Input:
%   info_bits : K-bit row vector (e.g., [1 0 1 1])
% Output:
%   x         : Encoded N-bit polar codeword

    % -------------------------------
    % Step 1: Reliability sequence (for N ≤ 32)
    % -------------------------------
    Q = [1 2 3 5 9 17 33 4 6 65 10 7 18 11 19 8 12 13 14 15 16 20 21 22 23 24 25 26 27 28 29 30 31 32];
    
    % -------------------------------
    % Step 2: Determine K and N
    % -------------------------------
    K = length(info_bits);
    
    % Ensure N > K (even when K is a power of 2)
    if bitand(K, K - 1) == 0
        N = 2^(log2(K) + 1);
    else
        N = 2^nextpow2(K);
    end
    n = log2(N);
    
    % Reliability sequence for this N
    Q1 = Q(Q <= N);
    
    % -------------------------------
    % Step 3: Insert frozen bits
    % -------------------------------
    u = zeros(1, N);
    info_pos = Q1(N - K + 1 : end);  % most reliable bit positions
    u(info_pos) = info_bits;          % fill info bits
    
    fprintf('Initial u (with frozen bits):\n');
    disp(u);
    
    % -------------------------------
    % Step 4: Butterfly structure encoding
    % -------------------------------
    m = 1; % block size
    for stage = 1:n
        for i = 1:2*m:N
            a = u(i : i + m - 1);           % upper branch
            b = u(i + m : i + 2*m - 1);     % lower branch
            
            % Perform polar transformation: [a ⊕ b, b]
            u(i : i + 2*m - 1) = [mod(a + b, 2), b];
        end
        fprintf('\nAfter stage %d:\n', stage);
        disp(u);
        m = m * 2; % double block size each stage
    end
    
    % -------------------------------
    % Step 5: Encoded output
    % -------------------------------
    x = u;  % final encoded bits
    fprintf('\nFinal encoded codeword:\n');
    disp(x);

end

% -------------------------------
% Example usage
% -------------------------------
msg = [1 0 1 1 0 1 0 0];    % input info bits (you can change this)
x = polar_encode_butterfly(msg);
