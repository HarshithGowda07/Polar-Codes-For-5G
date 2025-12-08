function crc = crc24a(msg_bits)
% CRC-24A polynomial: 0x1864CFB (x^24 + x^23 + x^18 + x^17 + x^14 + x^11 + x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x + 1)
poly = hex2dec('1864CFB');
reg = uint32(0);
for i=1:length(msg_bits)
    bit = uint32(msg_bits(i));
    msb = bitxor(bit, bitget(reg,24));
    reg = bitshift(reg,1);
    if msb
        reg = bitxor(reg, poly);
    end
    reg = bitand(reg, uint32(2^24-1));
end
crc = zeros(1,24);
for i=1:24
    crc(i) = bitget(reg, 25-i);
end
end
