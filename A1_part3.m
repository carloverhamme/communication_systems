%% Part 3b
clear; clc;

% Question 1: non-systematic generator matrix of rank 3
% (no set of three columns forms the identity: the column [1 0 0]' never appears)
G = [ 1  0  0  1  1  0 ; 1  0  1  0  1  1 ; 0  1  0  1  1  1 ];

k = size(G, 1); % = 3
n = size(G, 2); % = 6

assert(gfrank(G, 2) == k, 'G does not have rank 3 over F2');
fprintf('G has rank %d over F2\n\n', gfrank(G, 2));

% Question 2: all 2^k = 8 codewords
info_words = de2bi(0:(2^k - 1), k, 'left-msb'); % information vectors
C = mod(info_words * G, 2);                     % codewords
fprintf('Codewords:\n');
for idx = 1:size(C, 1)
    fprintf('   c_%d = [%s]  (v = [%s], weight %d)\n', idx-1, num2str(C(idx, :)), ...
        num2str(info_words(idx, :)), sum(C(idx, :)));
end
fprintf('d_min = %d\n\n', min(sum(C(2:end, :), 2)));

% Question 3: all 2^n = 64 received vectors
all_y = de2bi(0:(2^n - 1), n, 'left-msb');

% Question 4: decision regions with minimum-distance decoding
% (in case of a tie, y is included in all the corresponding regions)
D = cell(size(C, 1), 1);
n_ties = 0;                                     % received vectors with a tie
for i = 1:size(all_y, 1)
    y = all_y(i, :);
    dist = sum(xor(repmat(y, size(C, 1), 1), C), 2);
    min_dist = min(dist);
    closest_idx = find(dist == min_dist);
    if numel(closest_idx) > 1
        n_ties = n_ties + 1;
    end
    for idx = closest_idx'
        D{idx} = [D{idx}; y];
    end
end

% Question 5: write all the decision regions
for idx = 1:size(C, 1)
    c = C(idx, :);
    fprintf('=== Decision Region D(c_%d) for c = [%s] (Size: %d) ===\n', ...
        idx-1, num2str(c), size(D{idx}, 1));
    for r = 1:size(D{idx}, 1)
        fprintf('   [%s]  (distance %d)\n', num2str(D{idx}(r, :)), sum(xor(D{idx}(r, :), c)));
    end
    fprintf('\n');
end
fprintf('Received vectors with a tie: %d of %d\n\n', n_ties, size(all_y, 1));

% Question 6: verify D(c) = c + D(0) for all codewords
assert(all(C(1, :) == 0), 'The first codeword is not the all-zero codeword');
D0 = D{1};                                      % c_0 = [0 0 0 0 0 0]
all_valid = true;
for idx = 1:size(C, 1)
    c = C(idx, :);
    c_plus_D0 = mod(repmat(c, size(D0, 1), 1) + D0, 2);
    if ~isequal(sortrows(D{idx}), sortrows(c_plus_D0))
        all_valid = false;
        fprintf('Mismatch for codeword %d!\n', idx-1);
    end
end

if all_valid
    fprintf('D(c) = c + D(0) holds for all codewords!\n');
end

% Part 3c - Decision regions of a nonlinear code C(4,2)
 
% Question 1: four distinct codewords of length 4, including the all-zero one
C_nl = [0 0 0 0 ; 0 0 1 1 ; 0 1 0 1 ; 1 1 1 0];
M    = size(C_nl, 1);                           % number of codewords = 4
n_nl = size(C_nl, 2);                           % codeword length = 4
 
fprintf('\n\n===== Part 3c: nonlinear code =====\n');
fprintf('Codewords:\n');
for idx = 1:M
    fprintf('   c_%d = [%s]  (weight %d)\n', idx-1, num2str(C_nl(idx, :)), sum(C_nl(idx, :)));
end
 
% Check that the code is not linear: look for sums of codewords that are not codewords
is_linear = true;
for a = 1:M
    for b = a+1:M
        s = mod(C_nl(a, :) + C_nl(b, :), 2);
        if ~ismember(s, C_nl, 'rows')
            is_linear = false;
            fprintf('   c_%d + c_%d = [%s] is not a codeword\n', a-1, b-1, num2str(s));
        end
    end
end
if is_linear
    fprintf('The code is linear\n');
else
    fprintf('The code is NOT linear (not closed under addition)\n');
end
 
% Minimum distance (over all pairs of distinct codewords)
dmin_nl = inf;
for a = 1:M
    for b = a+1:M
        dmin_nl = min(dmin_nl, sum(xor(C_nl(a, :), C_nl(b, :))));
    end
end
fprintf('d_min = %d\n\n', dmin_nl);
 
% Question 2: decision regions with minimum-distance decoding (ties in all regions)
all_y4 = de2bi(0:(2^n_nl - 1), n_nl, 'left-msb');
D_nl = cell(M, 1);
n_ties_nl = 0;
for i = 1:size(all_y4, 1)
    y = all_y4(i, :);
    dist = sum(xor(repmat(y, M, 1), C_nl), 2);
    closest_idx = find(dist == min(dist));
    if numel(closest_idx) > 1
        n_ties_nl = n_ties_nl + 1;
    end
    for idx = closest_idx'
        D_nl{idx} = [D_nl{idx}; y];
    end
end
 
% Question 3: write the decision regions
for idx = 1:M
    c = C_nl(idx, :);
    fprintf('=== Decision Region D(c_%d) for c = [%s] (Size: %d) ===\n', ...
        idx-1, num2str(c), size(D_nl{idx}, 1));
    for r = 1:size(D_nl{idx}, 1)
        fprintf('   [%s]  (distance %d)\n', num2str(D_nl{idx}(r, :)), sum(xor(D_nl{idx}(r, :), c)));
    end
    fprintf('\n');
end
fprintf('Received vectors with a tie: %d of %d\n\n', n_ties_nl, size(all_y4, 1));
 
% Question 4: check whether D(c) = c + D(0) holds
D0_nl = D_nl{1};                                % c_0 = [0 0 0 0]
for idx = 1:M
    c = C_nl(idx, :);
    c_plus_D0 = mod(repmat(c, size(D0_nl, 1), 1) + D0_nl, 2);
    if isequal(sortrows(D_nl{idx}), sortrows(c_plus_D0))
        fprintf('c_%d: D(c) = c + D(0) HOLDS\n', idx-1);
    else
        fprintf('c_%d: D(c) = c + D(0) DOES NOT HOLD\n', idx-1);
        fprintf('      c + D(0)     = {%s}\n', strjoin(cellstr(num2str(sortrows(c_plus_D0), '%d')), ', '));
        fprintf('      D(c)         = {%s}\n', strjoin(cellstr(num2str(sortrows(D_nl{idx}), '%d')), ', '));
    end
end
 