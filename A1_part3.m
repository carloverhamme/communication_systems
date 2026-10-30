%% Part 3b
clear; clc;

% Question 1
G = [ 1  1  1  1  0  0 ; 0  1  1  0  1  0 ; 1  1  0  0  0  1 ];

k = size(G, 1); % =3
n = size(G, 2); % =6

% Question 2
info_words = de2bi(0:(2^k - 1), k, 'left-msb'); % infowords
C = mod(info_words * G, 2); % codewords
% display(C);

% Question 3
all_y = de2bi(0:(2^n - 1), n, 'left-msb'); % received vbectors

% Question 4
D = cell(size(C, 1), 1);
for i = 1:size(all_y, 1)
    y = all_y(i, :);
    dist = sum(xor(repmat(y, size(C, 1), 1), C), 2);
    min_dist = min(dist);
    closest_idx = find(dist == min_dist);
    for idx = closest_idx'
        D{idx} = [D{idx}; y];
    end
end

% Question 5
for idx = 1:size(C, 1)
    c = C(idx, :);
    fprintf('=== Decision Region D(c_%d) for c = [%s] (Size: %d) ===\n', ...
        idx-1, num2str(c), size(D{idx}, 1));
    for r = 1:size(D{idx}, 1)
        fprintf('   [%s]\n', num2str(D{idx}(r, :)));
    end
    fprintf('\n');
end

% Question 6
D0 = D{1}; % c_0 is [0 0 0 0 0 0]
all_valid = true;
for idx = 1:size(C, 1)
    c = C(idx, :);
    c_plus_D0 = mod(repmat(c, size(D0, 1), 1) + D0, 2);
    if ~isequal(sortrows(D{idx}), sortrows(c_plus_D0))
        all_valid = false;
        fprintf('Mismatch for codeword %d!', idx);
    end
end

if all_valid
    fprintf('D(c) = c + D(0) holds for all codewords!');
end

%---------------
% Part 3c