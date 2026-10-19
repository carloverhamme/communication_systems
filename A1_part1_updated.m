n=10;
k=5;

G1=[1,0,0,0,0, 1, 1, 1, 0, 0;0, 1, 0, 0, 0, 0, 1, 1, 1, 0; 0, 0, 1, 0, 0, 0, 0, 1, 1, 1;0,0,0, 1, 0,1, 0, 0, 1, 1;0, 0, 0, 0, 1,1, 1,0, 0, 1];
G2=[1,0,0,0,0, 1, 1, 0, 0, 0;0, 1, 0, 0, 0, 0, 1, 1, 1, 0; 0, 0, 1, 0, 0, 0, 0, 1, 1, 1;0,0,0, 1, 0,1, 0, 0, 1, 1;0, 0, 0, 0, 1,1, 1,0, 0, 1];

% PART 1
%Part 1.1
ac = 2^k; % Amount of codewords
inf_words = de2bi(0:ac-1, k, 'left-msb');

%Part 1.2
code_wordsG1 = mod(inf_words * G1, 2);
code_wordsG2 = mod(inf_words * G2, 2);

%Part 1.3
function weight = w_c(code_word) % w_c(codeword) computes hamming_weight
    weight = sum(code_word(:));
end


w_G1 = arrayfun(@(i) w_c(code_wordsG1(i, :)), 1:ac);
w_G2 = arrayfun(@(i) w_c(code_wordsG2(i, :)), 1:ac);

% Part 1.4
function Ai = multiplicity(codewords, i, ac)
    if i < 0 || i > 10
        Ai = 0;
        display("i should be in [0,10]");
        return
    end
    Ai = 0;
    for j = 1:ac
        if w_c(codewords(j, :)) == i
            Ai = Ai + 1;
        end
    end
end

A_G1 = arrayfun(@(i) multiplicity(code_wordsG1, i, ac), 0:n); % Multiplicity matrix of G1
A_G2 = arrayfun(@(i) multiplicity(code_wordsG2, i, ac), 0:n); % Multiplicity matrix of G2

% Part 1.5
d_min_G1 = min(w_G1(2:end));
d_min_G2 = min(w_G2(2:end));

A_min_G1 = A_G1(d_min_G1+1);
A_min_G2 = A_G2(d_min_G2+1);



% PART 2
% Part 2.1

% Part 2.1

function res = P_UE(p, A_i, n)
    res = 0;
    for i = 1:n
        res = res + A_i(i + 1) .* p.^i .* (1 - p).^(n - i);
    end
end

p = logspace(-6, -1);
P_UE_G1 = P_UE(p, A_G1, n);
P_UE_G2 = P_UE(p, A_G2, n);

% Part 2.2

loglog(p, P_UE_G1, p, P_UE_G2);
grid on;
legend('P_{UE, G1}', 'P_{UE, G2}');
xlabel('Probability p');
ylabel('Probability of an undetected error P_{UE}')
title('Ex. 1, Part 2 – Undetected error probability of codes G1 and G2');
% Part 2.3

lead_G1 = A_min_G1 .* p.^d_min_G1;
lead_G2 = A_min_G2 .* p.^d_min_G2;

figure;
loglog(p, P_UE_G1, 'b-', p, lead_G1, 'b*',p, P_UE_G2, 'r-', p, lead_G2, 'r*');
grid on;
legend('G1 exact', 'G1 leading term', 'G2 exact', 'G2 leading term');
xlabel('Probability p');
ylabel('Undetected error probability');
title('Ex. 1, Part 2 – Exact P_{UE} vs. leading term A_{min} p^{d_{min}}');
% Part 2.4

% G1 (blue) has better error-detection performance as the undetected error
% probability is lower across the whole p range.



% PART 3
% Part 3.1

% Part 3.1
H_T_G1 = gen2par(G1)';

% Part 3.2

y1 = code_wordsG1(10, :); % Existing codeword, random taken
error = [zeros(1, 9), 1]; % Error term
y2 = mod(y1 + error, 2); % Invalid codeword

s1 = mod(y1 * H_T_G1, 2);
s2 = mod(y2 * H_T_G1, 2);


% Part 3.3
U = [0 0 0 0 1
    0 1 0 0 0
    0 0 1 0 0
    0 0 0 1 0
    1 0 0 0 0]; % Swap first and 5th syndrome bits

H_T1 = mod(H_T_G1 * U, 2);
s1_H1 = mod(y1 * H_T1, 2);
s2_H1 = mod(y2 * H_T1, 2);


% Part 3.4
H_T2 = mod([H_T1, H_T1(:,1:2)], 2); % Zero syndrome condition stays the same here, just redundant columns are added

s1_H2 = mod(y1 * H_T2, 2);
s2_H2 = mod(y2 * H_T2, 2);


% PART 4
rng(1);
p_mc = [0.1, 0.09,0.08,0.07,0.06];
Pue10 = zeros(size(p_mc));
Pue100 = zeros(size(p_mc));
for i = 1:numel(p_mc)
    p_i = p_mc(i);
    ue = 0;
    transmitted = 0;

    while ue < 10
        rand_int = randi([1,32]);
        v = inf_words(rand_int, :);
        c = mod(v*G1,2);
        transmitted = transmitted+1;
        y = bsc(c, p_mc(i));
        s = mod(y * H_T_G1, 2);
        if all(s == 0) && any(y ~= c)
            ue = ue + 1;
        end
    end

    Pue10(i) = ue/transmitted;

    while ue < 100
        rand_int = randi([1,32]);
        v = inf_words(rand_int, :);
        c = mod(v*G1,2);
        transmitted = transmitted+1;
        y = bsc(c, p_mc(i));
        s = mod(y * H_T_G1, 2);
        if all(s == 0) && any(y ~= c)
            ue = ue + 1;
        end
    end

    Pue100(i) = ue/transmitted;
end

Pue_Formula = P_UE(p_mc, A_G1, n);

figure
loglog(p_mc, Pue_Formula, 'g-o', ...
     p_mc, Pue10, 'b-s', ...
     p_mc, Pue100, 'r-^')
grid on
xlabel('Probability p')
ylabel('P_{UE}')
legend('Analytical', 'Monte Carlo, 10 UEs', 'Monte Carlo, 100 UEs')
title('Ex. 1, Part 4 – Monte Carlo vs. analytical P_{UE} (code G1)');


% PART 5

p = logspace(-2, -1);

p0 = (1 - p).^n; % p0= P(COR, i=1)

p1 = zeros(size(p)); % p1 =P(UE, i=1)
for d = 1:n
    if A_G1(d+1) > 0
        p1 = p1 + A_G1(d+1) * (p.^d) .* ((1 - p).^(n - d));
    end
end

p2 = 1 - p0 - p1; % p2= P(REJ,i=1)




P_COR_1 = p0;
P_COR_2 = p0 .* (1 + p2);
P_COR_3 = p0 .* (1 + p2 + p2.^2);
P_err_1 = 1 - P_COR_1;
P_err_2 = 1 - P_COR_2;
P_err_3 = 1 - P_COR_3;

%Plot 1 - P(COR, i)
figure;
loglog(p, P_err_1, 'b-', 'DisplayName', 'i = 1');
hold on;
loglog(p, P_err_2, 'r-', 'DisplayName', 'i = 2');
loglog(p, P_err_3, 'g-',  'DisplayName', 'i = 3');
grid on;

xlabel('Probability p');
ylabel('1 - P(COR, i)');
legend;
title('Ex. 1, Part 5 – ARQ: 1 - P(COR, i) for i = 1, 2, 3 (code G1)');

%Part 6
EN = [ones(size(p)); 1 + p2; 1 + p2 + p2.^2; 1 ./ (1 - p2)];   % i = 1, 2, 3, inf


fprintf('G1: A_i = %s, d_min = %d, A_min = %d\n', mat2str(A_G1), d_min_G1, A_min_G1);
fprintf('G2: A_i = %s, d_min = %d, A_min = %d\n', mat2str(A_G2), d_min_G2, A_min_G2);
disp('H^T ='); disp(H_T_G1); disp('H1^T ='); disp(H_T1); disp('H2^T ='); disp(H_T2);
fprintf('y1 = %s | s(H) = %s, s(H1) = %s, s(H2) = %s\n', mat2str(y1), mat2str(s1), mat2str(s1_H1), mat2str(s1_H2));
fprintf('y2 = %s | s(H) = %s, s(H1) = %s, s(H2) = %s\n', mat2str(y2), mat2str(s2), mat2str(s2_H1), mat2str(s2_H2));
disp(table(p_mc', Pue_Formula', Pue10', Pue100', 'VariableNames', {'p','Analytical','MC10','MC100'}));
fprintf('p = %.2f: E[N] = %s\n', p(1), mat2str(EN(:,1)', 4)); fprintf('p = %.2f: E[N] = %s\n', p(end), mat2str(EN(:,end)', 4));