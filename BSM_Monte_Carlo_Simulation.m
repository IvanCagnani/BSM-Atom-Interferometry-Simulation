%% Monte Carlo simulation of the discovery of a BSM signal with an atom interferometer
%
% MATLAB code for the paper "A Background-Free Search for Physics Beyond the Standard Model Using Atom Interferometry" 
%
% Ivan Cagnani 
% 
% Linneaus University - Kalmar, Sweden
%
% ivan.cagnani@physics.org
% 
% November 16, 2025
%
% This script performs a Monte Carlo simulation to determine the
% sensitivity of the proposed experiment assuming a
% SUB-QUANTUM-LIMITED interferometer using spin squeezing.
%
% --- Methodology ---
% 1.  A 4-point differential quadrature is used to algebraically
%     isolate signals based on their parity (k-parity and V-parity),
%     matching the model in the accompanying paper.
% 2.  The interferometer noise is assumed to be 10x below the
%     Standard Quantum Limit (SQL) due to spin squeezing.
% 3.  The Monte Carlo simulation runs `N_trials` (e.g., 1,000,000)
%     single-shot experiments (`N_cycles_per_experiment = 1`) to find the
%     1-sigma noise floor for a *single 4-point measurement cycle* (`sigma_1`).
% 4.  This `sigma_1` is then analytically scaled by `1/sqrt(N_integration)`
%     to project the final sensitivity.




clear; close all; clc;




%% 1. Define Experimental Parameters
% --- Noise Model (Squeezed, Sub-SQL) ---
% Assumes a simultaneous dual-species interferometer with 20dB squeezing.
phase_noise_QPN = 0.0001; % (rad) 0.1 mrad (100 urad), 10x below SQL
% The noise on the simultaneous *difference* (Rb87-Rb85) is the
% quadrature sum of their individual QPNs.
noise_per_shot_simultaneous = sqrt(phase_noise_QPN^2 + phase_noise_QPN^2);
% noise_per_shot_simultaneous = sqrt(2) * 0.0001 rad ~= 0.141 mrad




% --- Input Signals (Amplitudes) ---
% These are the "true" physical signals we are simulating.
signal_BSM_amplitude = 1e-6;   % (Phi_BSM) 1.0 µrad (Our target signal)
systematic_V_amplitude = 1e-4;   % (Phi_V) 100 µrad (V-odd, k-even systematic)
systematic_k_amplitude = 5e-5;   % (Phi_k) 50 µrad (k-odd, V-even systematic)
% Note: Per the paper, Phi_0 (static) and Phi_kV (coupled) are also
% cancelled by the quadratures and are not simulated here.




%% 2. Define Monte Carlo Parameters
%
% We find the 1-sigma noise for *one* 4-point cycle, so we set
% N_cycles_per_experiment = 1.
%
% We run N_trials to build a smooth histogram of this single-shot noise.
% 1,000,000 trials gives a 0.1% precision (1/sqrt(N)) on the noise.
%
N_cycles_per_experiment = 1;
N_trials = 1000000; % Number of "virtual experiments" to run


% --- Define Target Integration Time ---
% This is used for the final sensitivity projection plot & analysis
N_integration_cycles = 126000;
cycle_time_sec = 60; % 60 seconds per 4-point cycle
total_time_days = (N_integration_cycles * cycle_time_sec) / (3600 * 24);




% Pre-allocate arrays for the result of each trial
quadrature_BSM_channel = zeros(N_trials, 1);
quadrature_k_channel   = zeros(N_trials, 1);




fprintf('Starting Monte Carlo Simulation (Squeezed Light Model)...\n');
fprintf('Running %d trials to find single-shot noise floor...\n', N_trials);




%% 3. Run Monte Carlo Simulation (Parallelized)


% --- Define Physical Model (from paper Sec. 2.4) ---
%
% The general model for the measured phase is (from Eq. 13):
% Phi_meas(k, V) = Phi_0 + k*(Phi_I + Phi_k) + V*(Phi_BSM + Phi_V) + (k*V)*Phi_kV
%
% The 4-point quadratures algebraically cancel Phi_0 and Phi_kV.
%
% This simulation focuses on noise and systematic isolation. We therefore
% omit the main g-dependent phase (Phi_I) as it is not needed to
% test the quadrature's ability to isolate the BSM and systematic
% channels from each other.
%
% We define the total amplitudes for the two non-cancelled channels:
% S_k_total (k-odd, V-even)  = Phi_k
% S_V_total (k-even, V-odd)  = Phi_BSM + Phi_V
%
% The simulation correctly measures these two sums in their
% respective quadratures.
%
% --- END MODEL DEFINITION ---


% Define the total amplitudes for the k-odd and V-odd channels
S_k_total = systematic_k_amplitude;
S_V_total = signal_BSM_amplitude + systematic_V_amplitude;


parfor i = 1:N_trials
    
    % --- Simulate one 4-point cycle (N_cycles_per_experiment = 1) ---
    
    % Generate the quantum-enhanced QPN noise for all 4 measurement types.
    noise_A = randn * noise_per_shot_simultaneous;
    noise_B = randn * noise_per_shot_simultaneous;
    noise_C = randn * noise_per_shot_simultaneous;
    noise_D = randn * noise_per_shot_simultaneous;




    % --- CALCULATE CORRECTED PHASES ---
    % Calculate the "measured" phase for each of the 4 shot types
    % (Signal + Systematics + Noise) based on the paper's model
    phase_A = (+1)*S_k_total + (+1)*S_V_total + noise_A; % A = Phi(+k, +V)
    phase_B = (-1)*S_k_total + (+1)*S_V_total + noise_B; % B = Phi(-k, +V)
    phase_C = (+1)*S_k_total + (-1)*S_V_total + noise_C; % C = Phi(+k, -V)
    phase_D = (-1)*S_k_total + (-1)*S_V_total + noise_D; % D = Phi(-k, -V)




    % --- Perform Quadrature Analysis (The Lock-In) ---
    % We combine the 4 measurements to isolate each component
    % based on the paper's quadratures (Sec 2.5 & 2.6).
    
    % BSM Channel (k-even, V-odd): (A+B) - (C+D) = 4*(S_BSM + S_V)
    % This is S_BSM from Eq. 24 in the paper.
    quadrature_BSM_channel(i) = (phase_A + phase_B) - (phase_C + phase_D);
    
    % k-Systematic Channel (k-odd, V-even): (A-B) + (C-D) = 4*S_k
    % This is S_I from Eq. 25 in the paper.
    quadrature_k_channel(i) = (phase_A - phase_B) + (phase_C - phase_D);
end




% Normalize all results by 4 to recover the original amplitudes
quadrature_BSM_channel = quadrature_BSM_channel / 4;
quadrature_k_channel   = quadrature_k_channel / 4;




fprintf('...Simulation complete.\n\n');




%% 4. Plot Results: Single-Shot Noise Distribution
% Use tiledlayout for better spacing and alignment
figure('Color', 'white', 'Position', [100 100 1000 500]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


% --- Plot the histogram of the *BSM channel* ---
nexttile(1); % Replaces subplot(1, 2, 1)


[N, edges] = histcounts(quadrature_BSM_channel * 1e6, 100, 'Normalization', 'pdf');
bar(edges(1:end-1) + diff(edges)/2, N, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');
hold on;
% Fit and plot a Gaussian
pd = fitdist((quadrature_BSM_channel * 1e6), 'Normal');
x_fit = linspace(min(edges), max(edges), 200);
y_fit = pdf(pd, x_fit);
p1 = plot(x_fit, y_fit, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Gaussian Fit');


% This channel measures (S_BSM + S_V)
true_signal_val_bsm = S_V_total * 1e6;
xline(true_signal_val_bsm, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off'); % Hide from legend


% Set title with formula as a subtitle
title_str_bsm = {
    'BSM Signal Channel (1-Shot)', ...
    sprintf('True Signal ($S_{BSM} + S_V$) = %.1f µrad', true_signal_val_bsm)
};
title(title_str_bsm);
xlabel('Measured Signal (µrad)');
ylabel('Probability Density');
grid on;
legend(p1, 'Location', 'northeast'); % Only show legend for the fit
xlim(true_signal_val_bsm + 5 * [-std(quadrature_BSM_channel), std(quadrature_BSM_channel)] * 1e6);




% --- Plot the histogram of the *k-Systematic channel* ---
nexttile(2); % Replaces subplot(1, 2, 2)


[N, edges] = histcounts(quadrature_k_channel * 1e6, 100, 'Normalization', 'pdf');
% Use a green color for the bars for better contrast
bar(edges(1:end-1) + diff(edges)/2, N, 'FaceColor', [0.4660 0.6740 0.1880], 'EdgeColor', 'none');
hold on;
pd = fitdist((quadrature_k_channel * 1e6), 'Normal');
x_fit = linspace(min(edges), max(edges), 200);
y_fit = pdf(pd, x_fit);
p2 = plot(x_fit, y_fit, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Gaussian Fit');


true_signal_val_k = S_k_total * 1e6;
xline(true_signal_val_k, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off'); % Hide from legend


% Set title with formula as a subtitle
title_str_k = {
    'k-Systematic Channel (1-Shot)', ...
    sprintf('True k-Systematic = %.1f µrad', true_signal_val_k)
};
title(title_str_k);
xlabel('Measured Signal (µrad)');
ylabel('Probability Density');
grid on;
legend(p2, 'Location', 'northeast'); % Only show legend for the fit
xlim(true_signal_val_k + 5 * [-std(quadrature_k_channel), std(quadrature_k_channel)] * 1e6);








%% 5. Plot Sensitivity vs. Integration Time
figure('Color', 'white', 'Position', [100 650 600 450]);


% Get 1-shot sigma from the simulation
sigma_1_shot_local = std(quadrature_BSM_channel);


% Create a log-spaced vector of integration cycles
integration_cycles_vec = logspace(0, log10(N_integration_cycles * 1.5), 100);


% Calculate the 1-sigma sensitivity (noise floor) vs. cycles
sensitivity_vec_1sigma = (sigma_1_shot_local ./ sqrt(integration_cycles_vec)) * 1e9; % in nrad


% Plot the sensitivity curve
loglog(integration_cycles_vec, sensitivity_vec_1sigma, 'b-', 'LineWidth', 2.5, 'DisplayName', '1-sigma Sensitivity (Noise Floor)');
hold on;


% Plot the target signal level
yline(signal_BSM_amplitude * 1e9, 'r--', 'LineWidth', 2, 'DisplayName', 'Target BSM Signal (1000 nrad)');


% Plot the 5-sigma discovery point
xline(N_integration_cycles, 'k:', 'LineWidth', 2, 'DisplayName', sprintf('5-sigma Detection Point (%.1f days)', (N_integration_cycles * 60) / (3600 * 24)));


title('Projected Sensitivity vs. Integration Time');
xlabel('Number of Integration Cycles');
ylabel('Projected Sensitivity (nrad)');
grid on;
legend('show', 'Location', 'southwest');
ylim([sensitivity_vec_1sigma(end)*0.1, sensitivity_vec_1sigma(1)*10]); % Set Y-limits
xlim([1, N_integration_cycles * 1.5]); % Set X-limits




%% 6. Final Analysis and Sensitivity Projection
% Get the 1-sigma noise floor for a single 4-point cycle
% The noise is identical in both channels, so we can use either.
sigma_1_shot = std(quadrature_BSM_channel);




fprintf('--- Monte Carlo Results (N = %d trials) ---\n', N_trials);
fprintf('Measured 1-Shot Noise Floor (sigma_1): %.3f µrad\n', sigma_1_shot * 1e6);




% --- Analytic Noise Floor Calculation (for comparison) ---
%
% 1. Noise per shot (dual-isotope differential):
%    noise_per_shot_simultaneous = sqrt(QPN_87^2 + QPN_85^2)
%                                = sqrt( (100 urad)^2 + (100 urad)^2 )
%                                = 141.4 urad
%
% 2. Error propagation for the final signal S_BSM = (A+B-C-D)/4:
%    Var(S_BSM) = (1/16) * (Var(A) + Var(B) + Var(C) + Var(D))
%               = (1/16) * (4 * Var(A))
%               = Var(A) / 4
%
% 3. Standard deviation (sigma_1) is sqrt(Var(S_BSM)):
%    sigma_1 = sqrt(Var(A) / 4) = sigma(A) / 2
%            = 141.4 urad / 2
%            = 70.7 urad
%
% This analytic value should match the simulated value.
analytic_sigma_1 = 0.5 * noise_per_shot_simultaneous;
fprintf('Analytic 1-Shot Noise Floor:         %.3f µrad\n', analytic_sigma_1 * 1e6);
fprintf(' (Note: These two numbers should be very close)\n\n');




% --- Sensitivity Projection ---
% Now we extrapolate to a long integration run.
% We use the *simulated* noise floor (sigma_1_shot) for the calculation.




% Final noise floor after N cycles of averaging
sigma_N_final = sigma_1_shot / sqrt(N_integration_cycles);




fprintf('--- Final Sensitivity Projection (Squeezed) ---\n');
fprintf('Target Integration:     %d cycles\n', N_integration_cycles);
fprintf('Equivalent Run Time:    %.1f days\n\n', total_time_days);
fprintf('Final 1-sigma Noise Floor (sigma_N): %.3e rad (%.3f nrad)\n', ...
        sigma_N_final, sigma_N_final * 1e9);




% Calculate the 1, 3, and 5-sigma sensitivity
sensitivity_1_sigma = 1 * sigma_N_final;
sensitivity_3_sigma = 3 * sigma_N_final;
sensitivity_5_sigma = 5 * sigma_N_final;




fprintf('1-sigma Sensitivity (Target): %.3f nrad\n', sensitivity_1_sigma * 1e9);
fprintf('3-sigma Sensitivity (Confidence): %.3f nrad\n', sensitivity_3_sigma * 1e9);
fprintf('5-sigma Sensitivity (Discovery):  %.3f nrad\n\n', sensitivity_5_sigma * 1e9);




fprintf('--- Sanity Check ---\n');
fprintf('Target BSM Signal:           %.1f nrad\n', signal_BSM_amplitude * 1e9);
fprintf('5-sigma Discovery Threshold: %.1f nrad\n', sensitivity_5_sigma * 1e9);
fprintf('Conclusion: The target signal is at the 5-sigma discovery level.\n');
