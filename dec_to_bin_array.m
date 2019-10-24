function bin_array = dec_to_bin_array(x,N)
% Convert a number to a binary array of 1s and 0s, ordered in ascending 
% bit. N specifies the number of places in the array (default is
% largest non-zero bit).
% Nate Zuk (2019)

if nargin<2,
    % Identify the maximum number of bits
    N = ceil(log2(x));
end

% Do the conversion (via long division)
bin_array = zeros(N,1);
for n = N:-1:1
    bit = 2^(n-1); % get the binary value at bit n
    % divide by that bit value
    bin_array(n) = floor(x/bit);
    % get the remainder
    x = x - bin_array(n)*bit;
end