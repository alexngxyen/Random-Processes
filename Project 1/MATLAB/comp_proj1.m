% Computer Project #1, Problem 1
% Date: February 9th, 2021
% Author: Alex Nguyen
% Description: Generate n = 100 i.i.d. Bernoulli RVs with P = 0.3

clc; clear; close all;

%% Initialize:
% Experiment Parameters
N = 1e4;     % Number of Bernoulli RV Sets (i.e., Experiments) 
n = 100;     % Number of Bernoulli RVs
p = 0.3;     % Probability of Success
k = 15:45;

% Theoretical Gaussian Distribution
y = -3:0.1:3;
my = 0;
sy = 1;
f_y = exp(-(y - my).^2./(2*sy^2))./(sy*sqrt(2*pi));

% Gaussian Approximation for Binomial RVs
P_k = @(k, n, p) 1/sqrt(2*pi*n*p*(1 - p))*exp(-(k - n*p).^2/(2*n*p*(1 - p)));

%% Bernoulli RV Distribution PMF:
% Generate Random Variables
b = rand(n, 1) < p;             % Approximate

% Plot
figure;
histogram(b, 'normalization', 'pdf');
xlabel('k'); ylabel('pmf');
title('Bernoulli Random Variable Distribution');

%% Sum of Bernoulli RVs (n = 100):
% Generate Bernoulli RVs
X = rand(N, n) < p;

% Sum of N Experiments
Sn = sum(X, 2);
Sn_true = P_k(k, n, p);

% Plot 
figure;
hold on;
histogram(Sn, 'normalization', 'pdf');
plot(k, Sn_true, 'r', 'linewidth', 2);
ylabel('pdf'); xlabel('k');
legend('Experimental', 'Theoretical');
title(sprintf('n = %d Bernoulli RV Sum', n)); 

%% Central Limit Theorem:
% Mean & Standard Deviation
mu = p;
sig2 = sqrt(p*(1 - p));

% Zero-Mean, Unit Variance Gaussian RV
Zn = (Sn - n*mu)./(sig2*sqrt(n));
Zn(Zn == 0) = [];

% Plot Histogram
figure;
hold on;
histogram(Zn, 'BinWidth', 0.2182179, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
title('Central Limit Theorem: Approximate Gaussian Distribution $\sim \mathcal{N}(0, 1)$', 'interpreter', 'latex');
legend('Approximate', 'Theoretical', 'fontsize', 10)
ylabel('pdf'); xlabel('x');

%% Central Limit Theorem For Different n's:
N1 = 5;
N2 = 50;

% Generate Bernoulli RVs
X1 = rand(N1, n) < p;
X2 = rand(N2, n) < p;

% Sum of N Experiments
SN1 = sum(X1, 2);
SN2 = sum(X2, 2);

% Zero-Mean, Unit Variance Gaussian RV
Zn1 = (SN1 - n*mu)./(sig2*sqrt(n));
Zn2 = (SN2 - n*mu)./(sig2*sqrt(n));

% Bin Width
w = @(n) 1/(sqrt(n*p*(1-p)));

% Plot 
figure;
subplot(3, 1, 1)
hold on;
histogram(Zn1, 'BinWidth', 1.2, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
ylabel('pdf'); xlabel('x');
legend(sprintf('n = %d', N1), 'Theoretical', 'location', 'best');
subplot(3, 1, 2)
hold on;
histogram(Zn2, 'BinWidth', .75, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
ylabel('pdf'); xlabel('x');
legend(sprintf('n = %d', N2), 'Theoretical', 'location', 'best');
subplot(3, 1, 3)
hold on;
histogram(Zn, 'BinWidth', 0.2182179, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
legend(sprintf('n = %d', n), 'Theoretical', 'location', 'best');
sgtitle('Central Limit Theorem: Approximate Gaussian Distribution $\sim \mathcal{N}(0, 1)$', 'interpreter', 'latex');
ylabel('pdf'); xlabel('x');
