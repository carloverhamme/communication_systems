%% Exercise 4a - Hamming (7,4) code: perfect code and analytical P_W(E)
clear; clc;

G = [1 0 0 0 1 1 0 ;
     0 1 0 0 1 0 1 ;
     0 0 1 0 0 1 1 ;
     0 0 0 1 1 1 1];
[k, n] = size(G);                               % k = 4, n = 7

% Codewords, minimum distance and error-correction capability
V = de2bi(0:(2^k - 1), k, 'left-msb');
C = mod(V * G, 2);
w = sum(C, 2);
dmin = min(w(2:end));                           % linear code: minimum nonzero weight
t = floor((dmin - 1) / 2);
fprintf('n = %d, k = %d, d_min = %d, t = %d\n', n, k, dmin, t);
fprintf('Weights of the codewords: %s\n\n', mat2str(sort(w')));

% Perfect code condition: 2^k * sum_{i=0}^{t} nchoosek(n, i) = 2^n
sphere_size = sum(arrayfun(@(i) nchoosek(n, i), 0:t));
fprintf('Vectors in a sphere of radius t: %d\n', sphere_size);
fprintf('2^k * sphere size = %d,  2^n = %d\n', 2^k * sphere_size, 2^n);
if 2^k * sphere_size == 2^n
    fprintf('-> The code is PERFECT\n\n');
else
    fprintf('-> The code is NOT perfect\n\n');
end

% Additional check: decision region D(0) obtained by exhaustive
% minimum-distance decoding must coincide with the sphere S(0, t)
Y = de2bi(0:(2^n - 1), n, 'left-msb');
D0 = [];
n_ties = 0;
for i = 1:size(Y, 1)
    dist = sum(xor(repmat(Y(i, :), size(C, 1), 1), C), 2);
    closest = find(dist == min(dist));
    if numel(closest) > 1
        n_ties = n_ties + 1;
    end
    if any(closest == 1)                        % first codeword is the all-zero one
        D0 = [D0; Y(i, :)];                     %#ok<AGROW>
    end
end
fprintf('|D(0)| = %d, max weight in D(0) = %d, ties = %d\n', ...
    size(D0, 1), max(sum(D0, 2)), n_ties);
if size(D0, 1) == sphere_size && max(sum(D0, 2)) <= t
    fprintf('-> D(0) = S(0, t)\n\n');
end

% Analytical word error probability for a perfect code
% P_W(E) = 1 - sum_{i=0}^{t} nchoosek(n, i) p^i (1-p)^(n-i)
p = logspace(-4, log10(0.5), 200);
P_correct = zeros(size(p));
for i = 0:t
    P_correct = P_correct + nchoosek(n, i) * p.^i .* (1 - p).^(n - i);
end
PW_hamming = 1 - P_correct;
PW_lead = nchoosek(n, t + 1) * p.^(t + 1);      % leading term for small p: 21 p^2

figure;
loglog(p, PW_hamming, 'b-', 'LineWidth', 1.5); hold on;
loglog(p, PW_lead, 'k--', 'LineWidth', 1);
grid on;
ylim([1e-7 1]);
xlabel('Crossover probability p');
ylabel('Word error probability P_W(E)');
legend('Analytical P_W(E)', 'Leading term 21p^2', 'Location', 'southeast');
title('Ex. 4a – Word error probability of the Hamming (7,4) code');
set(gca, 'FontSize', 12);

%% Exercise 4b - Monte Carlo simulation of P_W(E) for the Hamming (7,4) code
rng(1);                                         % reproducible results

n_wrong = 100;                                  % wrong codewords for each point
p_sim = [0.3 0.2 0.1 0.05 0.02 0.01 0.005 0.002];
[PW_sim, n_tx] = simulate_PW(G, C, p_sim, n_wrong);

% Analytical values at the simulated points
PW_th = zeros(size(p_sim));
for i = 0:t
    PW_th = PW_th + nchoosek(n, i) * p_sim.^i .* (1 - p_sim).^(n - i);
end
PW_th = 1 - PW_th;

fprintf('\nMonte Carlo simulation, Hamming (7,4), %d wrong codewords per point\n', n_wrong);
disp(table(p_sim', PW_th', PW_sim', n_tx', 100 * (PW_sim' - PW_th') ./ PW_th', ...
    'VariableNames', {'p', 'Analytical', 'Simulated', 'Transmitted', 'RelDiff_percent'}));

% Approximate 95% confidence intervals (same approach as Exercise 1, Part 4)
ci = 1.96 / sqrt(n_wrong) * PW_sim;

figure;
loglog(p, PW_hamming, 'b-', 'LineWidth', 1.5); hold on;
errorbar(p_sim, PW_sim, ci, 'ro', 'LineWidth', 1.2, 'MarkerSize', 6);
grid on;
ylim([1e-7 1]);
xlabel('Crossover probability p');
ylabel('Word error probability P_W(E)');
legend('Analytical', 'Monte Carlo (100 wrong codewords, 95% CI)', 'Location', 'southeast');
title('Ex. 4b – Analytical vs. simulated P_W(E), Hamming (7,4) code');
set(gca, 'FontSize', 12);

%% Exercise 4c - Extended Hamming (8,4) code
% Systematic generator matrix: Hamming (7,4) with an overall parity bit
G8 = [1 0 0 0 1 1 0 1 ;
      0 1 0 0 1 0 1 1 ;
      0 0 1 0 0 1 1 1 ;
      0 0 0 1 1 1 1 0];
[k8, n8] = size(G8);                            % k = 4, n = 8

V8 = de2bi(0:(2^k8 - 1), k8, 'left-msb');
C8 = mod(V8 * G8, 2);
w8 = sum(C8, 2);
dmin8 = min(w8(2:end));
t8 = floor((dmin8 - 1) / 2);
fprintf('\n\n===== Exercise 4c: (8,4) code =====\n');
fprintf('n = %d, k = %d, d_min = %d, t = %d\n', n8, k8, dmin8, t8);
fprintf('Weights of the codewords: %s\n', mat2str(sort(w8')));

% Perfect code condition
sphere8 = sum(arrayfun(@(i) nchoosek(n8, i), 0:t8));
fprintf('2^k * sphere size = %d,  2^n = %d\n', 2^k8 * sphere8, 2^n8);
if 2^k8 * sphere8 == 2^n8
    fprintf('-> The code is PERFECT\n\n');
else
    fprintf('-> The code is NOT perfect\n\n');
end

% Decision region D(0) by exhaustive minimum-distance decoding.
% With random tie-breaking, an error vector e in a tie among m codewords
% (including 0) is corrected with probability 1/m.
Y8 = de2bi(0:(2^n8 - 1), n8, 'left-msb');
a_w = zeros(1, n8 + 1);                         % a_w(w+1) = sum of 1/m over e in D(0) of weight w
count_w = zeros(1, n8 + 1);                     % vectors of D(0) for each weight
ties_w = zeros(1, n8 + 1);                      % tied vectors of D(0) for each weight
n_ties8 = 0;
for i = 1:size(Y8, 1)
    e = Y8(i, :);
    dist = sum(xor(repmat(e, size(C8, 1), 1), C8), 2);
    closest = find(dist == min(dist));
    m = numel(closest);
    if m > 1
        n_ties8 = n_ties8 + 1;
    end
    if any(closest == 1)                        % 0 is among the closest codewords
        we = sum(e);
        a_w(we + 1) = a_w(we + 1) + 1 / m;
        count_w(we + 1) = count_w(we + 1) + 1;
        ties_w(we + 1) = ties_w(we + 1) + (m > 1);
    end
end
fprintf('Received vectors with a tie: %d of %d\n', n_ties8, 2^n8);
fprintf('D(0) by weight:\n');
for wv = 0:n8
    if count_w(wv + 1) > 0
        fprintf('   weight %d: %2d vectors (%2d ties), probability of correct decoding sum = %g\n', ...
            wv, count_w(wv + 1), ties_w(wv + 1), a_w(wv + 1));
    end
end

% Exact P_W(E) with random tie-breaking, and Lecture 7 upper bound
PW8_exact = ones(size(p));
for wv = 0:n8
    PW8_exact = PW8_exact - a_w(wv + 1) * p.^wv .* (1 - p).^(n8 - wv);
end
PW8_bound = ones(size(p));
for i = 0:t8
    PW8_bound = PW8_bound - nchoosek(n8, i) * p.^i .* (1 - p).^(n8 - i);
end
fprintf('Max |P_W(8,4) - P_W(7,4)| over p: %g\n\n', max(abs(PW8_exact - PW_hamming)));

% Monte Carlo simulation (same points and stopping criterion as Part 4b)
rng(2);
[PW8_sim, n_tx8] = simulate_PW(G8, C8, p_sim, n_wrong);
PW8_th = ones(size(p_sim));                     % exact values at the simulated points
for wv = 0:n8
    PW8_th = PW8_th - a_w(wv + 1) * p_sim.^wv .* (1 - p_sim).^(n8 - wv);
end
PW8_ub = ones(size(p_sim));
for i = 0:t8
    PW8_ub = PW8_ub - nchoosek(n8, i) * p_sim.^i .* (1 - p_sim).^(n8 - i);
end
fprintf('Monte Carlo simulation, (8,4) code, %d wrong codewords per point\n', n_wrong);
disp(table(p_sim', PW8_th', PW8_ub', PW8_sim', n_tx8', 100 * (PW8_sim' - PW8_th') ./ PW8_th', ...
    'VariableNames', {'p', 'Exact', 'UpperBound', 'Simulated', 'Transmitted', 'RelDiff_percent'}));

% Comparison of the two codes
ci8 = 1.96 / sqrt(n_wrong) * PW8_sim;
figure;
loglog(p, PW_hamming, 'b-', 'LineWidth', 1.5); hold on;
loglog(p, PW8_exact, 'r--', 'LineWidth', 1.5);
loglog(p, PW8_bound, 'k:', 'LineWidth', 1.5);
errorbar(p_sim, PW_sim, ci, 'bo', 'LineWidth', 1.2, 'MarkerSize', 6);
errorbar(p_sim, PW8_sim, ci8, 'rs', 'LineWidth', 1.2, 'MarkerSize', 7);
grid on;
ylim([1e-7 1]);
xlabel('Crossover probability p');
ylabel('Word error probability P_W(E)');
legend('(7,4) analytical', '(8,4) analytical (random ties)', '(8,4) upper bound', ...
    '(7,4) Monte Carlo', '(8,4) Monte Carlo', 'Location', 'southeast');
title('Ex. 4c – Hamming (7,4) vs. extended Hamming (8,4)');
set(gca, 'FontSize', 12);

%% Local functions
function [PW_sim, n_tx] = simulate_PW(G, C, p_vec, n_wrong)
% Monte Carlo estimate of the word error probability over a BSC with
% minimum-distance decoding. For each p, random codewords are transmitted
% until n_wrong wrong decoded codewords are observed.
    [k, n] = size(G);
    batch = 20000;                              % codewords processed together
    PW_sim = zeros(size(p_vec));
    n_tx = zeros(size(p_vec));
    for ip = 1:numel(p_vec)
        p_i = p_vec(ip);
        errors = 0;
        transmitted = 0;
        while errors < n_wrong
            V = randi([0 1], batch, k);         % random information vectors
            X = mod(V * G, 2);                  % transmitted codewords
            E = rand(batch, n) < p_i;           % BSC error vectors
            Y = mod(X + E, 2);                  % received vectors
            X_hat = decode_md(Y, C);            % minimum-distance decoding
            wrong = find(any(X_hat ~= X, 2));   % indices of wrong codewords
            if errors + numel(wrong) >= n_wrong
                % stop exactly at the n_wrong-th wrong codeword
                last = wrong(n_wrong - errors);
                transmitted = transmitted + last;
                errors = n_wrong;
            else
                transmitted = transmitted + batch;
                errors = errors + numel(wrong);
            end
        end
        PW_sim(ip) = errors / transmitted;
        n_tx(ip) = transmitted;
    end
end

function X_hat = decode_md(Y, C)
% Minimum Hamming distance decoding of each row of Y.
% Distances: number of positions where y and c differ.
% Ties are broken at random (a small random number is added to each distance).
    dist = Y * (1 - C)' + (1 - Y) * C';         % size: (number of vectors) x (number of codewords)
    dist = dist + 0.1 * rand(size(dist));       % random tie-breaking
    [~, idx] = min(dist, [], 2);
    X_hat = C(idx, :);
end