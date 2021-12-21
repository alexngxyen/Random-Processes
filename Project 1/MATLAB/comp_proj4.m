% Computer Project #1, Problem 4
% Date: February 9th, 2021
% Author: Alex Nguyen
% Description: Verify law of large numbers by calculating sample means.

clc; clear; close all;

%% Initialize Number of Random Experiments:
n_low = 5;     % Low Number 
n_mid = 100;   % Medium Number 
n_high = 1e4;  % High Number 

%% Bernoulli RV:
% Initialize
p = 0.4;       % Probability of Success

% Low Number of Trials
X1 = rand(n_low, 1) < p;
mu_1 = mean(X1);

% Medium Number of Trials
X2 = rand(n_mid, 1) < p;
mu_2 = mean(X2);

% Large Number of Trials
X3 = rand(n_high, 1) < p;
mu_3 = mean(X3);

% Print Results
fprintf(['Bernoulli Random Variable with Expected Value, E[X] = %4.2f: \n\t' ...
    '   i) Low Number of Trials: Mn = %4.3f (N = %d)\n\t  ii) Medium Number of Trials: ' ...
    'Mn = %4.3f (N = %d)\n\t iii) High Number of Trials: Mn = %4.3f (N = %d)' ...
    '\n\n'], p, mu_1, n_low, mu_2, n_mid, mu_3, n_high);

%% Binomial RV:
N = 6;         % Bernoulli Trials
p = 0.3;       % Probability of Success

% Low Number of Trials
X1 = sum(rand(n_low, N) < p, 2);
mu_1 = mean(X1);

% Large Number of Trials
X2 = sum(rand(n_mid, N) < p, 2);
mu_2 = mean(X2);

% Large Number of Trials
X3 = sum(rand(n_high, N) < p, 2);
mu_3 = mean(X3);

% Print Results
fprintf(['Binomial Random Variable with Expected Value, E[X] = %4.2f: \n\t' ...
    '   i) Low Number of Trials: Mn = %4.3f (N = %d) \n\t  ii) Medium Number of Trials '...
    ': Mn = %4.3f (N = %d) \n\t iii) High Number of Trials: Mn = %4.3f (N = %d)' ...
    '\n\n'], N*p, mu_1, n_low, mu_2, n_mid, mu_3, n_high);

%% Law of Large Numbers:
fprintf(['The law of large numbers (LLN) states the sample mean (Mn) will converge' ...
    ' in probability to the expected value (or theoretical mean) as\nthe number of' ...
    ' samples in a sequence increases. Above, the two examples illustrate this phenomena.' ...
    ' Thus, the LLN is verified!\n\n']);