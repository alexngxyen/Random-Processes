% Computer Project #1, Problem 5
% Date: February 9th, 2021
% Author: Alex Nguyen
% Description: Estimate mean of X^2 where X ~ N(0, 1)

clc; clear; close all;

%% Initialize:
% Different Number of Gaussian RV Realizations
n1 = 1e1;                 % Low Number of Trials        
n2 = 5e3;                 % Medium Number of Trials  
n3 = 1e5;                 % High Number of Trials 

% Gaussian RVs
X1 = randn(n1, 1);        
X2 = randn(n2, 1);     
X3 = randn(n3, 1);      

%% var(X) = E[X^2] - (E[X]).^2  --> E[X^2] = var(X) + (E[X]).^2 :
% Sample Mean Estimate
mu_hat1 = sum(X1)/n1;
mu_hat2 = sum(X2)/n2;
mu_hat3 = sum(X3)/n3;

% Unbiased Variance Estimate
sig2_hat1 = sum((X1 - mu_hat1).^2)/(n1 - 1);
sig2_hat2 = sum((X2 - mu_hat2).^2)/(n2 - 1);
sig2_hat3 = sum((X3 - mu_hat3).^2)/(n3 - 1);

% E[X^2] Estimate
EX2_hat1 = sig2_hat1 + mu_hat1.^2;
EX2_hat2 = sig2_hat2 + mu_hat2.^2;
EX2_hat3 = sig2_hat3 + mu_hat3.^2;

%% Estimate for E[X^2]:
fprintf(['Estimated Expected value of X^2 (or E[X^2]) with mu = 1: \n\t '...
    '  i) Low Number of Samples: E[X^2] = %4.4f (N = %d) \n\t  ii) Medium Number of Samples: E[X^2] = %4.4f (N = %d) \n\t' ...
    ' iii) High Number of Samples: E[X^2] = %4.4f (N = %d) \n\n'], EX2_hat1, n1, EX2_hat2, n2, EX2_hat3, n3);

