clear all;
close all;
clc;

%% Load Dataset
data = load('doa_dataset_correlated.mat');
obserX = data.X;
N = size(obserX, 1);                  % Number of array elements
L = size(obserX, 2);                  % Number of snapshots
M = 2;                                % Number of sources
phoIndex = 1;                         % Selected correlation index
trialIndex = 1;                       % Selected trial

% Extract observation matrix for selected trial and rho
obserVector = obserX(:, :, trialIndex, phoIndex);

% Compute Sample Covariance Matrix (uncentered or centered)
covarianceMatrix = (obserVector * obserVector') / L;

%% Run DOA Algorithms
thetaDegree = -90 : 0.1 : 90;

% 1. MUSIC Algorithm
[theta_hat_music, musicSpectrum] = music_doa(covarianceMatrix, N, M, thetaDegree, data.dOverLambda);

% 2. Maximum Likelihood (ML) Algorithm
theta_hat_ml = ml_doa(covarianceMatrix, N, M, thetaDegree, data.dOverLambda);

%% Convert MUSIC Spectrum to Decibels (dB) and Normalize
musicSpectrum_dB = 10 * log10(musicSpectrum);
musicSpectrum_dB = musicSpectrum_dB - max(musicSpectrum_dB); % Normalize peak to 0 dB

%% Enhanced Plotting
figure('Color', 'w', 'Position', [100, 100, 900, 550]);

% Plot Normalized MUSIC Pseudospectrum
p1 = plot(thetaDegree, musicSpectrum_dB, 'k-', 'LineWidth', 1.8, 'DisplayName', 'MUSIC Pseudospectrum');
hold on; grid on; grid minor;

% Plot True DOAs
trueDOAs = data.thetaTrue;
for k = 1:length(trueDOAs)
    xline(trueDOAs(k), 'r--', sprintf('True: %.1f^o', trueDOAs(k)), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', ...
        'LabelHorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.8 0 0], ...
        'HandleVisibility', iff(k == 1, 'on', 'off'), 'DisplayName', 'True DOAs');
end

% Plot MUSIC Estimated Peaks
music_peaks_dB = zeros(1, length(theta_hat_music));
for k = 1:length(theta_hat_music)
    [~, idx] = min(abs(thetaDegree - theta_hat_music(k)));
    music_peaks_dB(k) = musicSpectrum_dB(idx);
end
p2 = plot(theta_hat_music, music_peaks_dB, 'ro', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'MUSIC Detected Peaks');

% Plot ML Estimates as Vertical Stem/Lines
for k = 1:length(theta_hat_ml)
    xline(theta_hat_ml(k), 'm-.', sprintf('ML: %.1f^o', theta_hat_ml(k)), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.7 0 0.7], ...
        'HandleVisibility', iff(k == 1, 'on', 'off'), 'DisplayName', 'ML Estimates');
end

% Peak Text Annotations
for k = 1:length(theta_hat_music)
    text(theta_hat_music(k), music_peaks_dB(k) + 1.5, ...
        sprintf('MUSIC: %.1f^o', theta_hat_music(k)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9, 'Color', 'k');
end

% Graph Details
title(sprintf('DOA Estimation: MUSIC vs. Maximum Likelihood (M = %d, \\rho = %.2f)', ...
    N, data.rho(phoIndex)), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Angle of Arrival \theta (Degrees)', 'FontSize', 11);
ylabel('Normalized MUSIC Spectrum (dB)', 'FontSize', 11);
xlim([-90 90]);
ylim([min(musicSpectrum_dB)-5, 5]);
legend('Location', 'southwest', 'FontSize', 10);

%% Print Results to Command Window
fprintf('\n======================================================\n');
fprintf('DOA Estimation Results (rho = %f):\n', data.rho(phoIndex));
fprintf('======================================================\n');
fprintf('True DOAs (deg):        [%s ]\n', sprintf(' %.2f', trueDOAs));
fprintf('MUSIC Estimates (deg):  [%s ]\n', sprintf(' %.2f', theta_hat_music));
fprintf('ML Estimates (deg):     [%s ]\n', sprintf(' %.2f', theta_hat_ml));
fprintf('======================================================\n');


%% ====================================================================
%% Helper & DOA Functions
%% ====================================================================

function [theta_hat, musicSpectrum] = music_doa(Rhat, N, M, theta_grid, dOverLambda)
    % Eigenvalue Decomposition of Covariance Matrix
    [V, D] = eig(Rhat);
    eigenValues = diag(D);
    
    [~, sortIndex] = sort(eigenValues, 'descend');
    Vsorted = V(:, sortIndex);
    
    % Extract Noise Subspace (N - M dimensions)
    noiseSubspace = Vsorted(:, M+1:end);
    
    theta = theta_grid * (pi / 180);
    musicSpectrum = zeros(1, length(theta));
    
    % Compute MUSIC Pseudospectrum
    for k = 1:length(theta)
        steeringVector = exp(1j * 2 * pi * dOverLambda * sin(theta(k)) * (0:N-1)');
        deno = steeringVector' * (noiseSubspace * noiseSubspace') * steeringVector;
        musicSpectrum(k) = 1 / abs(deno);
    end
    
    % Peak Detection
    [pks, locs] = findpeaks(musicSpectrum);
    [~, peakOrder] = sort(pks, 'descend');
    numSources = min(M, numel(peakOrder));
    
    theta_hat = NaN(1, M);
    theta_hat(1:numSources) = sort(theta_grid(locs(peakOrder(1:numSources))), 'ascend');
end

function theta_hat = ml_doa(Rhat, N, M, grid, dOverLambda)
    num_grid_points = numel(grid);
    A_all = zeros(N, num_grid_points);
    grid_rad = grid * (pi / 180);
    
    for k = 1:num_grid_points
        A_all(:, k) = exp(+1j * 2 * pi * dOverLambda * sin(grid_rad(k)) * (0:N-1)');
    end
    
    max_cost = -Inf;
    coarse_winner_idx = [1, 2];
    
    for i = 1:num_grid_points
        for j = (i + 1):num_grid_points
            a_i = A_all(:, i);
            a_j = A_all(:, j);
            A_cand = [a_i, a_j];
            
            if rcond(A_cand' * A_cand) < 1e-4
                continue;
            end
            
            B = A_cand' * Rhat * A_cand;
            cost = real(trace((A_cand' * A_cand) \ B));
            
            if cost > max_cost
                max_cost = cost;
                coarse_winner_idx = [i, j];
            end
        end
    end
    
    theta_coarse = grid(coarse_winner_idx);
    refine_range1 = (theta_coarse(1) - 1):0.1:(theta_coarse(1) + 1);
    refine_range2 = (theta_coarse(2) - 1):0.1:(theta_coarse(2) + 1);
    
    max_cost_fine = -Inf;
    theta_hat = theta_coarse;
    
    for k1 = 1:numel(refine_range1)
        th1 = refine_range1(k1);
        a_fine1 = exp(+1j * 2 * pi * dOverLambda * sin(th1 * pi / 180) * (0:N-1)');
        
        for k2 = 1:numel(refine_range2)
            th2 = refine_range2(k2);
            if abs(th1 - th2) < 0.1, continue; end
            
            a_fine2 = exp(+1j * 2 * pi * dOverLambda * sin(th2 * pi / 180) * (0:N-1)');
            A_fine = [a_fine1, a_fine2];
            
            if rcond(A_fine' * A_fine) < 1e-4, continue; end
            
            B_fine = A_fine' * Rhat * A_fine;
            cost_fine = real(trace((A_fine' * A_fine) \ B_fine));
            
            if cost_fine > max_cost_fine
                max_cost_fine = cost_fine;
                theta_hat = [th1, th2];
            end
        end
    end
    
    theta_hat = sort(theta_hat, 'ascend');
end

% Inline Conditional Helper Function
function val = iff(cond, trueVal, falseVal)
    if cond, val = trueVal; else, val = falseVal; end
end