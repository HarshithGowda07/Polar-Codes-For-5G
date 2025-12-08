function [crc_len, crc_func_handle] = polar_crc_select(K)
% Select CRC length and function based on payload K (5G-like rules)
% returns length (6/11/24) and a function handle to compute CRC bits

    if K <= 19
        crc_len = 6;
        crc_func_handle = @crc6;
    elseif K <= 105
        crc_len = 11;
        crc_func_handle = @crc11;
    else
        crc_len = 24;
        crc_func_handle = @crc24a;
    end
end
