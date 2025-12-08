% polar_simulation_main_dynamic.m
% Main script: enter payload (msg bits) OR set random_payload = true.
% Automatically chooses CRC (6/11/24), N (next pow2 >= K+CRC), default L, runs SC and CA-SCL, plots BER.

clc; clear; close all;

%% ===== USER OPTIONS =====
% Provide payload bits as a row vector OR set random_payload = true
random_payload = false;   % if true, payload will be randomized each frame for BER sim
% Example fixed message (uncomment and edit if you want fixed payload):
% payload = [1 0 1 1 0 1 0 1 0 1 1 1 0 0 1 0 0 0 1 0 1 0 1 1 1 0];

% Choose simulation mode:
doBERsim = true;    % true: run BER vs SNR (Monte Carlo), false: single-shot test
numFrames = 200;    % frames per SNR (reduce for debugging)
SNRdBlist = -2:2:8; % Eb/N0 list for BER plot

% Default list size (you can override later)
default_L = 8;

% Max N allowed (you set 1024)
N_max = 1024;

%% ===== Get payload length / message =====
if ~exist('payload','var') || isempty(payload)
    if ~random_payload
        % user didn't set payload; ask to type payload here
        % For automated run, fallback to a short random payload:
        payload = randi([0 1],1,120); % fallback 120 bits (you can edit)
    else
        payload = []; % will randomize each frame
    end
end
K_msg = length(payload);

fprintf('User payload length K_msg = %d\n', K_msg);

%% ===== Choose CRC length dynamically (6/11/24) =====
[crc_len, crc_func_handle] = polar_crc_select(K_msg);
fprintf('Selected CRC length = %d bits\n', crc_len);

K_total = K_msg + crc_len;

%% ===== Choose N = smallest power of 2 >= K_total, but <= N_max =====
N = 2^nextpow2(K_total);
if N > N_max
    error('Required N = %d exceeds N_max = %d. Increase N_max or reduce payload.', N, N_max);
end
n = log2(N);
fprintf('Chosen blocklength N = %d (n=%d), total K(with CRC) = %d\n', N, n, K_total);

%% ===== Choose reliability sequence Q (5G NR truncated) =====
Q = get_5g_polar_sequence(N);  % returns vector 1..N in reliability order
info_pos = Q(N-K_total+1:end);
frozen_pos = Q(1:N-K_total);

%% ===== Choose list size L based on K_total (simple heuristic) =====
if K_total <= 128
    L = 8;
elseif K_total <= 256
    L = 16;
else
    L = 32;
end
% allow override by default_L
if exist('default_L','var') && ~isempty(default_L)
    L = default_L;
end
fprintf('Using list size L = %d\n', L);

%% ===== Run either single-shot test or BER simulation =====
if ~doBERsim
    % Single-shot encode/decode
    if isempty(payload)
        payload = randi([0 1], 1, K_msg);
    end
    crc_bits = crc_func_handle(payload);
    info_bits = [payload crc_bits]; % length K_total
    x = polar_encoder(info_bits, Q);
    % Channel
    EbNo_dB = 4; Rate = K_total / N;
    sigma2 = 1/(2*Rate*10^(EbNo_dB/10));
    r = (1-2*x) + sqrt(sigma2)*randn(1,N);
    Lx = 2*r/sigma2;
    % SC decode
    [u_sc, msg_sc] = polar_sc_decode_x_then_u(Lx, Q, K_total);
    % CA-SCL decode (Tal-Vardy)
    [u_list, msg_cascl] = polar_scl_talvardy(Lx, Q, K_total, crc_len, crc_func_handle, L);
    fprintf('payload original:\n'); disp(payload);
    fprintf('SC-decoded payload:\n'); disp(msg_sc(1:K_msg));
    fprintf('CA-SCL-decoded payload:\n'); disp(msg_cascl(1:K_msg));
    return;
end

% BER simulation
BER_SC = zeros(size(SNRdBlist));
BER_CASCL = zeros(size(SNRdBlist));

for sidx = 1:length(SNRdBlist)
    EbNo_dB = SNRdBlist(sidx);
    EbNo = 10^(EbNo_dB/10);
    Rate = K_total / N;
    sigma2 = 1/(2*Rate*EbNo);
    sigma = sqrt(sigma2);

    bit_errors_sc = 0;
    bit_errors_cascl = 0;
    total_bits = 0;

    for frame = 1:numFrames
        % payload: either fixed or random each frame
        if random_payload
            payload_f = randi([0 1], 1, K_msg);
        else
            payload_f = payload;
        end
        crc_bits = crc_func_handle(payload_f);
        info_bits = [payload_f crc_bits]; % length K_total
        x = polar_encoder(info_bits, Q);
        s = 1 - 2*x;
        r = s + sigma*randn(1,N);
        Lx = 2*r / sigma2;

        % SC decode
        [u_sc, msg_sc] = polar_sc_decode_x_then_u(Lx, Q, K_total);
        payload_hat_sc = msg_sc(1:K_msg);
        bit_errors_sc = bit_errors_sc + sum(payload_f ~= payload_hat_sc);

        % CA-SCL decode
        [u_cascl, msg_cascl] = polar_scl_talvardy(Lx, Q, K_total, crc_len, crc_func_handle, L);
        payload_hat_cascl = msg_cascl(1:K_msg);
        bit_errors_cascl = bit_errors_cascl + sum(payload_f ~= payload_hat_cascl);

        total_bits = total_bits + K_msg;
    end

    BER_SC(sidx) = bit_errors_sc / total_bits;
    BER_CASCL(sidx) = bit_errors_cascl / total_bits;
    fprintf('SNR=%2d dB | SC BER=%.3g | CA-SCL(L=%d) BER=%.3g\n', EbNo_dB, BER_SC(sidx), L, BER_CASCL(sidx));
end

% Plot
figure; semilogy(SNRdBlist, BER_SC, '-o','LineWidth',2); hold on;
semilogy(SNRdBlist, BER_CASCL, '-s','LineWidth',2); grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER'); title(sprintf('SC vs CA-SCL (N=%d,K=%d,crc=%d,L=%d)',N,K_msg,crc_len,L));
legend('SC','CA-SCL');
