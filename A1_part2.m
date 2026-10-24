clear; clc;

% Point 1: primitive polynomial p(D) = D^4 + D + 1, degree r = 4
r = 4;
p = [1 0 0 1 1];                   

% Point 2: generator polynomial g(D) = (D + 1) p(D)
g = mod(conv([1 1], p), 2);         % product of polynomials over F2
crc_poly = 'z^5 + z^4 + z^2 + 1';   % the same g(D), in MATLAB syntax
fprintf('p(D) = %s\n', poly_str(p));
fprintf('g(D) = %s  (degree %d)\n\n', poly_str(g), numel(g) - 1);

% Points 3-8: start from k = 2 and increase k until weight-2 codewords appear
k = 2;
k_list = [];                        % values of k considered
A_list = {};                        % weight distribution for each k
dmin_list = [];                     % minimum distance for each k

while true
    n  = k + r + 1;                 % codeword length
    nc = 2^k;                       % number of codewords

    % All information vectors (one per row) and systematic CRC encoding
    V = de2bi(0:nc-1, k, 'left-msb');
    encoder = comm.CRCGenerator('Polynomial', crc_poly);
    C = zeros(nc, n);
    for m = 1:nc
        C(m, :) = encoder(V(m, :)')';
    end

    % Sanity checks: systematic encoding, every codeword multiple of g(D)
    assert(isequal(C(:, 1:k), V), 'Encoding is not systematic');
    for m = 1:nc
        [~, rem_poly] = deconv(C(m, :), g);
        assert(all(mod(rem_poly, 2) == 0), 'Codeword not divisible by g(D)');
    end

    % Hamming weights and weight distribution: A(i+1) = A_i, i = 0..n
    w = sum(C, 2);
    A = accumarray(w + 1, 1, [n + 1, 1])';
    dmin = min(w(2:end));           % first row is the all-zero codeword

    % Point 5: table for k = 2
    if k == 2
        fprintf('Table for k = 2 (n = %d)\n', n);
        fprintf('%-6s %-10s %-26s %-10s %s\n', 'v', 'v(D)', 'c(D)', 'c', 'w_H(c)');
        for m = 1:nc
            fprintf('%-6s %-10s %-26s %-10s %d\n', sprintf('%d', V(m, :)), ...
                poly_str(V(m, :)), poly_str(C(m, :)), sprintf('%d', C(m, :)), w(m));
        end
        fprintf('\n');
    end

    k_list(end + 1) = k;            %#ok<SAGROW>
    A_list{end + 1} = A;            %#ok<SAGROW>
    dmin_list(end + 1) = dmin;      %#ok<SAGROW>

    % Point 8: stop when codewords of weight 2 appear
    if A(3) > 0
        fprintf('Weight-2 codewords appear for k = %d (n = %d):\n', k, n);
        for m = find(w == 2)'
            fprintf('  c = %s,  c(D) = %s\n', sprintf('%d', C(m, :)), poly_str(C(m, :)));
        end
        fprintf('\n');
        break
    end
    k = k + 1;
end

% Points 6 and 9: weight distribution A_i for each value of k
n_max = k_list(end) + r + 1;
fprintf('Weight distribution A_i\n');
fprintf('%3s %3s %5s |', 'k', 'n', 'dmin');
fprintf('%5d', 0:n_max);
fprintf('\n');
for t = 1:numel(k_list)
    A = A_list{t};
    fprintf('%3d %3d %5d |', k_list(t), k_list(t) + r + 1, dmin_list(t));
    fprintf('%5d', [A, zeros(1, n_max + 1 - numel(A))]);
    fprintf('\n');
end

%% Local functions
function s = poly_str(bits)
% Polynomial in D from a coefficient vector (highest degree first)
    deg = numel(bits) - 1;
    terms = {};
    for idx = 1:numel(bits)
        if bits(idx) == 1
            e = deg - idx + 1;
            if e == 0
                terms{end + 1} = '1';               %#ok<AGROW>
            elseif e == 1
                terms{end + 1} = 'D';               %#ok<AGROW>
            else
                terms{end + 1} = sprintf('D^%d', e); %#ok<AGROW>
            end
        end
    end
    if isempty(terms)
        s = '0';
    else
        s = strjoin(terms, ' + ');
    end
end
%% 

%-----------
% PART 2b


crc_enc = comm.CRCGenerator('Polynomial', crc_poly);
crc_dec = comm.CRCDetector('Polynomial', crc_poly);

% 2. Select an information vector (k = 2)
% Example: v = [1; 0]
v = [1; 0];
c = crc_enc(v);  % Systematic codeword of length n = 2 + 5 = 7

fprintf('--- Codeword ---\n');
fprintf('Transmitted Codeword c: [%s]\n', num2str(c'));

% 3. Define Error Patterns (column vectors of length 7)
% e1: Burst of length L = 5 (indices 1 to 5) -> Guaranteed detected
e1 = [1; 0; 0; 0; 1; 0; 0];

% e2: Burst equal to g(D) of length L = 6 -> Undetected
e2 = [1; 1; 0; 1; 0; 1; 0];

% 4. Received vectors
y1 = mod(c + e1, 2);
y2 = mod(c + e2, 2);

% 5. CRC Syndrome Detection
[~, err1] = crc_dec(y1);
[~, err2] = crc_dec(y2);

fprintf('\n--- Detection Results ---\n');
fprintf('e1 (Burst length 5): err = %d -> %s\n', err1, ternary(err1 == 1, 'DETECTED (Syndrome != 0)', 'UNDETECTED'));
fprintf('e2 (Burst length 6): err = %d -> %s\n', err2, ternary(err2 == 1, 'DETECTED', 'UNDETECTED (Syndrome == 0)'));

function out = ternary(cond, valTrue, valFalse)
if cond, out = valTrue; else, out = valFalse; end
end