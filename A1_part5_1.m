% Part 5_1
clear; clc;

n = 63; k = 57; r = n - k; % r = 6

% Ex.1
all_cols = de2bi(1:(2^r - 1), r, 'left-msb')';
is_unit = (sum(all_cols, 1) == 1);
P_T = all_cols(:, ~is_unit);
H = [P_T, eye(r)];
G = [eye(k), P_T'];

weights_vec = 2.^(r-1:-1:0)';
syn_dec_cols = weights_vec' * H;
lut = zeros(1, 2^r - 1);
for col = 1:n
    lut(syn_dec_cols(col)) = col;
end

% Ex 2
p_analytical = logspace(-4, -0.7, 300);
Pw_analytical = 1 - ((1 - p_analytical).^n + n .* p_analytical .* (1 - p_analytical).^(n - 1));

% EX. 4
p_sim = [1e-3, 2e-3, 4e-3, 7e-3, 1e-2, 1.5e-2, 2e-2];
Pw_sim = zeros(size(p_sim));
target_errors = 100;

for idx = 1:length(p_sim)
    p = p_sim(idx);
    errors = 0;
    frames = 0;
    while errors < target_errors
        frames = frames + 1;
        v = randi([0 1], 1, k);
        c = mod(v * G, 2);
        noise = (rand(1, n) < p);
        y = mod(c + noise, 2);

        % Bounded-distance syndrome decoding (t = 1)
        s = mod(y * H', 2);
        s_val = s * weights_vec;
        c_hat = y;
        if s_val > 0
            c_hat(lut(s_val)) = mod(c_hat(lut(s_val)) + 1, 2);
        end

        if ~isequal(c_hat, c)
            errors = errors + 1;
        end
    end
    Pw_sim(idx) = errors / frames;
end

% Ex. 5 and 6 (also 3)
figure;
loglog(p_analytical, Pw_analytical, 'b-', 'LineWidth', 1.8); hold on;
loglog(p_sim, Pw_sim, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);
grid on;
xlabel('Channel Bit Error Probability p');
ylabel('Codeword Error Probability P_W(E)');
title('Hamming (63,57) Bounded-Distance Syndrome Decoding over BSC');
legend('Analytical P_W(E)', 'Monte Carlo Simulation (100 errors)', 'Location', 'northwest');

