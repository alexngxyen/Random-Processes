% Computer Project #1, Problem 3
% Date: February 9th, 2021
% Author: Alex Nguyen
% Description: Generate CDF of a Gaussian RV with N(2, 3) from a uniform RV
% in [0, 1].

clc; clear; close all;

%% Intialize:
% Experiment Parameters
n = 5e3;                 % Uniform RV Sample Realizations
mu = 2;                  % Mean
sig2 = 3;                % Variance

% Theoretical Gaussian Distribution
y = -4:0.1:8;
my = mu;
sy = sqrt(sig2);
f_y = exp(-(y - my).^2./(2*sy^2))./(sy*sqrt(2*pi));

%% Generate Gaussian Random Variables:
% Two Independent Uniform RVs
Z1 = rand(n, 1); 
Z2 = rand(n, 1);

% Zero-Mean, Unit-Variance X ~ N(0, 1)
R2 = -2*log(Z1);                       % Exponential RV [exp(-R^2/2)]
TH = 2*pi*Z2;                          % Uniformly Distributed RV [0, 2*pi]
Z = sqrt(R2).*sin(TH);
X = sqrt(sig2)*Z + mu;                 % Gaussian RV

%% Plot:
% Gaussian (Normal) Distribution PDF [X ~ N(2, 3)]
figure; 
hold on;
histogram(X, 'Normalization', 'pdf');
plot(y, f_y, 'r', 'linewidth', 2)
hold off;
xlabel('x'); ylabel('pdf');
legend('Aproximate', 'Theoretical', 'location', 'best');
title('Gaussian (Normal) Distribution PDF [$\sim \mathcal{N}(2, 3)$]', 'interpreter', 'latex');

% Gaussian (Normal) Distribution CDF [X ~ N(2, 3)]
figure;
ecdf(X);
title('Gaussian (Normal) Distribution Empirical CDF [$\sim \mathcal{N}(2, 3)$]', 'interpreter', 'latex');
