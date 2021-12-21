% Computer Project #2, Problem 1
% Date: March 10th, 2021
% Author: Alex Nguyen
% Description: Simulate autoregressive (AR) random process.

clc; clear; close all;

%% Initialize:
Y0 = 0;           % Initial Y Value
L = 500;          % Number of Discrete Values in Time Index
t = 1:L;          % Time Index Values (disc-time)
N = 3;            % Number of Realizations

%% Autoregressive (AR) Random Process, Yn = alpha*Yn-1 + Xn:
% Variance Function
sig_x = @(x) 1 - x.^2;

% Auto-Regressive Variance Parameters
alp1 = 0.3;   sig_x1 = sig_x(alp1);  % Process 1
alp2 = 0.95;  sig_x2 = sig_x(alp2);  % Process 2

% % Textbook (Figure 10.9) AR Variance Parameters
% alp1 = 0.1;   sig_x1 = sig_x(alp1);  % Process 1
% alp2 = 0.75;  sig_x2 = sig_x(alp2);  % Process 2

%% Sample Realizations:
% For, \alpha = 0.3
Yn1 = zeros(N, L);     % Preallocation
Y0_1 = Y0;

% For, \alpha = 0.95
Yn2 = zeros(N, L);     % Preallocation
Y0_2 = Y0;

for n = 1:N
    for ii = 1:L
        % Zero-Mean White Noise
        Xn1 = sqrt(sig_x1)*randn;
        Xn2 = sqrt(sig_x2)*randn;
        
        % Auto-Regressive Process
        Yn1(n, ii) = alp1*Y0_1 + Xn1;
        Yn2(n, ii) = alp2*Y0_2 + Xn2;
        
        % Update
        Y0_1 = Yn1(n, ii);
        Y0_2 = Yn2(n, ii);
    end
end

% Plot
fig = figure;
subplot(2, 1, 1)
hold on;
for ii = 1:N
    plot(t, Yn1(ii, :), 'linewidth', 2);
end
hold off;
legend('$\alpha = 0.3$', 'interpreter', 'latex');
xlim([t(1) t(end)]);

subplot(2, 1, 2)
hold on;
for ii = 1:N
    plot(t, Yn2(ii, :), 'linewidth', 2);
end
hold off;
legend('$\alpha = 0.95$', 'interpreter', 'latex');
xlim([t(1) t(end)]);

han = axes(fig, 'visible', 'off'); 
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, '$Y_n(n, \xi)$ ', 'interpreter', 'latex', 'fontsize', 14);
xlabel(han, 'n (# of Samples)', 'fontsize', 14);
sgtitle('Sample Realizations for the Auto-Regressive Model', 'fontweight', 'normal', 'fontsize', 12);
grid on;
    
%% Auto-Correlation, Ry(k):
% Estimated Auto-Correlation 
[Ry1_est, lag1] = xcorr(Yn1(1, :), 'biased');  % \alpha = 0.3
[Ry2_est, lag2] = xcorr(Yn2(1, :), 'biased');  % \alpha = 0.95

% Plot Estimated Ry(k)
figure;  
hold on;
plot(lag1, Ry1_est, 'r', 'linewidth', 2); 
plot(lag2, Ry2_est, 'b', 'linewidth', 2);
hold off;
xlabel('Lags', 'fontsize', 14);
ylabel('$R_Y(k)$', 'interpreter', 'latex', 'fontsize', 14);
legend('$\alpha = 0.3$', '$\alpha = 0.95$', 'interpreter', 'latex');
title('Estimated Auto-Correlation for Auto-Regressive Process', 'fontweight', 'normal', 'fontsize', 12);
grid on;

%% Power Spectral Density, Sy(k):
% Initialize
n = length(Ry1_est);     % Length of Auto-Correlation Vector
f = [-n/2:n/2 - 1]./n;   % Normalized Frequency [Hz]

% Estimated Zero-Centered Power Spectral Densities
Sy1_fft_est = fft(Ry1_est);            
Sy1_Pavg = Sy1_fft_est.*conj(Sy1_fft_est) ./ sqrt(n);
Sy1_est = fftshift(Sy1_Pavg);

Sy2_fft_est = fft(Ry2_est);
Sy2_Pavg = Sy2_fft_est.*conj(Sy2_fft_est) ./ sqrt(n);
Sy2_est = fftshift(Sy2_Pavg);

% Plot with Output [dB/Hz]
figure;
hold on;
plot(f, 10*log10(Sy1_est), 'r', 'linewidth', 1.25);
plot(f, 10*log10(Sy2_est), 'b', 'linewidth', 1.25);
hold off;
ylabel('$S_Y(f)$ [dB/Hz]', 'interpreter', 'latex', 'fontsize', 14);
xlabel('Normalized Frequency, f [Hz]', 'fontsize', 14);
legend('$\alpha = 0.3$', '$\alpha = 0.95$', 'interpreter', 'latex');
title('Estimated Power Spectral Density for Auto-Regressive Process', 'fontweight', 'normal', 'fontsize', 12);
grid on;
xlim([f(1) f(end)]);

% Plot with Output [Hz^-1]
fig = figure;
subplot(2, 1, 1);
plot(f, Sy1_est, 'r', 'linewidth', 1.25);
legend('$\alpha = 0.3$', 'interpreter', 'latex');
xlim([f(1) f(end)]);
grid on;

subplot(2, 1, 2);
plot(f, Sy2_est, 'b', 'linewidth', 1.25);
legend('$\alpha = 0.95$', 'interpreter', 'latex');
xlim([f(1) f(end)]);
grid on;

han = axes(fig, 'visible', 'off'); 
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, '$S_Y(f)$ [Hz$^{-1}$]', 'interpreter', 'latex', 'fontsize', 14);
xlabel(han, 'Normalized Frequency, f [Hz]', 'fontsize', 14);
sgtitle('Estimated Power Spectral Density for Auto-Regressive Process', 'fontweight', 'normal', 'fontsize', 12);
grid on;

