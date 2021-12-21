% Computer Project #1, Problem 2
% Date: February 9th, 2021
% Author: Alex Nguyen
% Description: Generate n = 100 i.i.d. Poisson RVs with lambda = 0.3

clc; clear; close all;

%% Initialize:
% Poisson RV Parameters
N = 1e4;                 % Number of Random Experiments
n = 100;                 % Number of Poisson RVs
lam = 0.3;               % Rate Parameter
k = 0:5;                 % Observations

% Theoretical Poisson RV
p_k = @(k, a) (a.^k./factorial(k)).*exp(-a);

% Theoretical Gaussian Distribution
y = -3:0.1:3;
my = 0;
sy = 1;
f_y = exp(-(y - my).^2./(2*sy^2))./(sy*sqrt(2*pi));

%% Poisson RV Distribution PMF:
% Generate Random Variables
p = poissrnd(lam*ones(n, 1));  % Approximate
p_true = p_k(k, lam);          % True 

% Plot
figure;
hold on;
histogram(p, 'BinWidth', 1, 'normalization', 'pdf');
plot(k, p_true, 'r', 'linewidth', 2);
hold off;
xlabel('k'); ylabel('pmf');
legend('Approximate', 'Theoretical', 'location', 'best');
title('Poisson Random Variable Distribution');

%% Sum of Poisson RVs (n = 100):
% Generate Sum of Poisson RV
Sn = zeros(N, 1);              % Preallocation
 
for ii = 1:N
    Sn(ii) = sum(poissrnd(lam*ones(n, 1)));
end

% Plot
figure;
histogram(Sn, 'normalization', 'pdf');
xlabel('k'); ylabel('pdf')
title(sprintf('n = %d Poisson RV Sum', n));

%% Central Limit Theorem:
% Mean & Standard Deviation
mu = lam;
sig2 = sqrt(lam);

% Zero-Mean, Unit Variance Gaussian RV
Zn = (Sn - n*mu)./(sig2*sqrt(n));
Zn(Zn == 0) = [];

% Plot Histogram
figure;
hold on;
histogram(Zn, 'BinWidth', 0.19, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
title('Central Limit Theorem: Approximate Gaussian Distribution $\sim \mathcal{N}(0, 1)$', 'interpreter', 'latex');
legend('Approximate', 'Theoretical', 'fontsize', 10)
ylabel('pdf'); xlabel('x');

%% Central Limit Theorem For Different n's:
n1 = 5;
n2 = 50;

% Generate Bernoulli RVs
Sn1 = zeros(N, 1);        % Preallocation
Sn2 = Sn1;

for ii = 1:N
    Sn1(ii) = sum(poissrnd(lam*ones(n1, 1)));
    Sn2(ii) = sum(poissrnd(lam*ones(n2, 1)));
end

% Zero-Mean, Unit Variance Gaussian RV
Zn1 = (Sn1 - n1*mu)./(sig2*sqrt(n1));
Zn2 = (Sn2 - n2*mu)./(sig2*sqrt(n2));

% Plot 
figure;
subplot(3, 1, 1)
hold on;
histogram(Zn1, 'Binwidth', 1, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
ylabel('pdf'); xlabel('x');
legend(sprintf('n = %d', n1), 'Theoretical', 'location', 'best');
subplot(3, 1, 2)
hold on;
histogram(Zn2, 'Binwidth', .65, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
ylabel('pdf'); xlabel('x');
legend(sprintf('n = %d', n2), 'Theoretical', 'location', 'best');
subplot(3, 1, 3)
hold on;
histogram(Zn, 'BinWidth', 0.2182179, 'normalization', 'pdf');
plot(y, f_y, 'r', 'LineWidth', 2);
hold off;
legend(sprintf('n = %d', n), 'Theoretical', 'location', 'best');
sgtitle('Central Limit Theorem: Approximate Gaussian Distribution $\sim \mathcal{N}(0, 1)$', 'interpreter', 'latex');
ylabel('pdf'); xlabel('x');
