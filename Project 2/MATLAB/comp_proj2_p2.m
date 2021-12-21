% Computer Project #2, Problem 2
% Date: March 10th, 2021
% Author: Alex Nguyen
% Description: Simulate random phase random process model for 100
% realizations.

clc; clear; close all;

%% Initialize:
L = 30;             % Number of Discrete Values in Time Index
t = 1:L;            % Time Index Values (disc-time)
N = 100;            % Number of Realizations 

%% Sinusoid with Random Phase, X_n = cos(w*n + Omega):
% Random Phase Sinusoid Function
X_n_Xi = @(n, th) cos(0.2*pi*n + th);

% Preallocate
Xn = zeros(N, L);
Omega = zeros(N, 1);

for ii = 1:N
    % Uniform Random Variable [-pi, pi]
    Omega(ii) = (pi - (-pi))*rand + (-pi);
    
    % Single Realization of Random Phase Sinusoid
    Xn(ii, :) = X_n_Xi(t, Omega(ii));
end

%% Plot:
% 100 Different Realizations of Xn vs n
figure;
hold on;
for ii = 1:N
    plot(t, Xn(ii, :), '.', 'linewidth', 2);
end
hold off;
xlabel('n (# of Samples)', 'fontsize', 12);
ylabel('$X_n(n, \xi)$ ', 'interpreter', 'latex', 'fontsize', 12);
title('Sinusoid with Random Phase', 'fontsize', 14, 'fontweight', 'normal');
xlim([0 t(end)]);

% Different Number of Realizations of Random Sinusoid
fig = figure;
subplot(3, 1, 1)
plot(t, Xn(75, :), '-', 'linewidth', 2);
xlim([t(1) t(end)]);
title('1 Realization', 'fontsize', 14, 'fontweight', 'normal');

subplot(3, 1, 2)
hold on;
for ii = 1:5
    plot(t, Xn(ii, :), '-', 'linewidth', 2);
end
hold off;
xlim([t(1) t(end)]);
title('5 Realizations', 'fontsize', 14, 'fontweight', 'normal');

subplot(3, 1, 3)
hold on;
for ii = 1:50
    plot(t, Xn(ii, :), '-', 'linewidth', 2);
end
hold off;
xlim([t(1) t(end)]);
title('50 Realizations', 'fontsize', 14, 'fontweight', 'normal');

han = axes(fig, 'visible', 'off'); 
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, '$X_n(n, \xi)$ ', 'interpreter', 'latex', 'fontsize', 14);
xlabel(han, 'n (# of Samples)', 'fontsize', 14);

