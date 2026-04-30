%% ==========================================================================
% BSM Atom Interferometry: SSFM Physics Engine & Monte Carlo Validation
% Author: Ivan Cagnani
% License: MIT
% Date: 2026-04-30
% Description: Ab initio 1D Split-Step Fourier Method simulation of a
%              dual-isotope atom interferometer searching for anomalous
%              gravitational-electromagnetic scalar couplings. Includes
%              spin-squeezed state modeling, mid-fringe covariance extraction,
%              and true-seed Monte Carlo validation.
% Requirements: MATLAB R2020a+, Statistics and Machine Learning Toolbox
% Note: Includes 120 true seeds from Random.org
% Paper: Cagnani, I. (2026). A Background-Free Search for Physics Beyond the 
% Standard Model Using Atom Interferometry. Zenodo. 
% Paper link: https://doi.org/10.5281/zenodo.18896187
% Contact: ivan.cagnani@physics.org
% ===========================================================================
clear; clc; close all;

%% ========================================================================
% PART 1: THE PHYSICS ENGINE (Runs Once via Fast SSFM)
% =========================================================================
fprintf('--- PART 1: AB INITIO PHYSICS SIMULATION ---\n');
hbar = 1.0545718e-34;
amu  = 1.660539e-27;
m85  = 84.9118 * amu;
m87  = 86.9092 * amu;

g = 9.81;                 
T_drift = 0.01;           
k_eff = 1e6;              
Omega = 2 * pi * 1e5;     

% Accelerations
a_BSM_base = 2.0e-7;  
a_BSM_85 = 48 * a_BSM_base; 
a_BSM_87 = 50 * a_BSM_base; 

F_Stark = m87 * 5.0e-6;   
a_Stark_87 = F_Stark / m87;
a_Stark_85 = F_Stark / m85; 

L = 10e-3;            
N_grid = 1048576; % 2^20 SSFM Grid

fprintf('Initializing Fast SSFM Grid for N = %d...\n', N_grid);
dx = L/N_grid;
x_vals = dx * (0:N_grid-1)' - L/2;
dp = 2*pi/L;
p_vals = hbar * dp * [0:N_grid/2-1, -N_grid/2:-1]'; 

N_steps = 50; 
dt = T_drift / N_steps;
U_T85 = exp(-1i * (p_vals.^2 / (2 * m85 * hbar)) * dt);
U_T87 = exp(-1i * (p_vals.^2 / (2 * m87 * hbar)) * dt);

sigma = 1.5e-5; 
psi_init_g = exp(-0.25 * (x_vals.^2) / sigma^2);
psi_init_g = psi_init_g / sqrt(sum(abs(psi_init_g).^2 * dx)); 
psi_init_e = zeros(N_grid, 1);

True_Phases = zeros(4, 2); % [Shot, Isotope]

density_g_halfway = zeros(N_grid, 1); 
density_e_halfway = zeros(N_grid, 1); 

fprintf('Running Fast MZ Sequence to find True Physical Phases...\n');
tic;
for shot = 1:4
    if shot == 1 || shot == 3
        k_dir = 1; unwind = exp(-1i * k_eff * x_vals);
    else
        k_dir = -1; unwind = exp( 1i * k_eff * x_vals);
    end
    
    if shot == 1 || shot == 2
        V_sign = 1; 
    else
        V_sign = -1; 
    end
    
    a_85 = g + V_sign * (a_BSM_85 + a_Stark_85);
    a_87 = g + V_sign * (a_BSM_87 + a_Stark_87);
    
    U_V85_half = exp(-1i * (m85 * a_85 / hbar) * x_vals * (dt/2));
    U_V87_half = exp(-1i * (m87 * a_87 / hbar) * x_vals * (dt/2));
    
    % 85Rb
    g85 = psi_init_g; e85 = psi_init_e;
    [g85, e85] = apply_pulse(g85, e85, pi/2, k_dir * k_eff, x_vals);
    for s = 1:N_steps, [g85, e85] = apply_drift(g85, e85, U_T85, U_V85_half); end
    if shot == 1
        density_g_halfway = abs(g85).^2;
        density_e_halfway = abs(e85).^2;
    end
    [g85, e85] = apply_pulse(g85, e85, pi, k_dir * k_eff, x_vals);
    for s = 1:N_steps, [g85, e85] = apply_drift(g85, e85, U_T85, U_V85_half); end
    True_Phases(shot, 1) = angle(sum( conj(g85) .* e85 .* unwind ));
    
    % 87Rb
    g87 = psi_init_g; e87 = psi_init_e;
    [g87, e87] = apply_pulse(g87, e87, pi/2, k_dir * k_eff, x_vals);
    for s = 1:N_steps, [g87, e87] = apply_drift(g87, e87, U_T87, U_V87_half); end
    [g87, e87] = apply_pulse(g87, e87, pi, k_dir * k_eff, x_vals);
    for s = 1:N_steps, [g87, e87] = apply_drift(g87, e87, U_T87, U_V87_half); end
    True_Phases(shot, 2) = angle(sum( conj(g87) .* e87 .* unwind ));
end
fprintf('Physics Engine completed in %.2f seconds.\n\n', toc);

expected_BSM   = k_eff * (a_BSM_87 - a_BSM_85) * T_drift^2;
expected_Stark = k_eff * (a_Stark_87 - a_Stark_85) * T_drift^2;
target_signal  = expected_BSM + expected_Stark;

%% ========================================================================
% PART 2: THE NOISE ENGINE & TRUE-SEED MONTE CARLO VALIDATION
% =========================================================================
fprintf('--- PART 2: STATISTICAL MONTE CARLO (TRUE-SEED VALIDATION) ---\n');

% 1. Load True Random Seeds from random.org (10 sets of 12)
true_seeds = [
636728303, 799724275, 391091111, 489826994, 801537481, 100170805, 432781465, 79051795, 51089656, 55331488, 984892439, 494643647
405623237, 340846039, 30324471, 672764240, 872209380, 495624087, 743727585, 788967642, 843701274, 68002538, 5827836, 691443865
58969388, 847030943, 946249564, 825956110, 6617372, 60008627, 313700243, 418578344, 578911832, 115277805, 40455194, 837542236
461351755, 610649127, 426227994, 582467505, 326580215, 135580126, 88569740, 740312478, 942062008, 714850122, 522924940, 517389687
246682198, 491957456, 615309908, 755690824, 933559671, 72069570, 101684244, 101477631, 26022104, 651215332, 183704858, 714588702
346947571, 933030461, 36877250, 977703150, 464185415, 691754126, 520253774, 94609297, 415028907, 698789240, 928105152, 944727468
335476563, 871114114, 694065327, 990073192, 687530416, 83992793, 292657959, 76216511, 283859513, 897926534, 827975617, 576678401
92527827, 701437980, 942087418, 594488841, 846342701, 98462482, 532077409, 384900696, 760261379, 854843978, 360187684, 322659800
408443572, 5301684, 261292149, 589470608, 220707115, 474855366, 254352581, 629569372, 401034482, 19576903, 830632946, 521880758
273519581, 522896319, 588181803, 822263696, 285298133, 960424464, 103855070, 131429715, 525810402, 106827257, 347440318, 760563953
];

[num_sims, num_runs] = size(true_seeds);
total_iterations = num_sims * num_runs;
sigma_values = zeros(total_iterations, 1);
counter = 1;

% Monte Carlo Constants
N_cycles = 25000;          
n_points = 1000;          
N_atoms  = 1e6;           
sigma_diff_drift = 2.0e-6; 

tau_ref = 3^(1/6) / N_atoms^(2/3); 
tau_star = (2 / (3 * N_atoms))^(1/6) * tau_ref; 
C_tau = cos(tau_star)^(N_atoms - 1); 
sigma_z_sq = 2^(-1/3) * N_atoms^(-4/3); 

demo_z85 = []; demo_z87 = []; demo_dphi = 0;

fprintf('Executing %d Validation Runs (each with %d MC cycles)...\n', total_iterations, N_cycles);
tic;

for sim = 1:num_sims
    for run = 1:num_runs
        
        % Initialize the PRNG with the True Seed
        current_seed = true_seeds(sim, run);
        rng(current_seed, 'twister');
        
        S_BSM_results = zeros(N_cycles, 1);
        S_I_results   = zeros(N_cycles, 1);
        
        for cycle = 1:N_cycles
            extracted_diff_phases = zeros(4, 1);
            drift_rate = randn() * sigma_diff_drift;
            
            for shot = 1:4
                Phi_85_true = True_Phases(shot, 1);
                Phi_87_true = True_Phases(shot, 2);
                
                shot_drift = drift_rate * shot; 
                true_delta_phi = (Phi_87_true - Phi_85_true) + shot_drift;
                
                phi_cn = linspace(0, 2*pi, n_points)' + rand() * 2*pi;
                mid_fringe_offset = pi/2;
                
                phi_87_meas = phi_cn + Phi_87_true + shot_drift + mid_fringe_offset; 
                phi_85_meas = phi_cn + Phi_85_true;
                
                z_87 = -C_tau * sin(phi_87_meas) + randn(n_points, 1) * sqrt(sigma_z_sq);
                z_85 = -C_tau * sin(phi_85_meas) + randn(n_points, 1) * sqrt(sigma_z_sq);
                
                % Capture demo data only on the very first loop so Part 3 doesn't break
                if sim == 1 && run == 1 && cycle == 1 && shot == 1
                    demo_z85 = z_85; demo_z87 = z_87; demo_dphi = true_delta_phi;
                end
                
                est_phase = extract_differential_phase(z_85, z_87, C_tau);
                extracted_diff_phases(shot) = est_phase;
            end
            
            Phi_A = extracted_diff_phases(1);
            Phi_B = extracted_diff_phases(2);
            Phi_C = extracted_diff_phases(3);
            Phi_D = extracted_diff_phases(4);
            
            S_BSM_results(cycle) = 0.25 * (Phi_A - Phi_B - Phi_C + Phi_D);
            S_I_results(cycle)   = 0.25 * (Phi_A - Phi_B + Phi_C - Phi_D);
        end
        
        mean_BSM = mean(S_BSM_results);
        std_BSM  = std(S_BSM_results);
        
        % Calculate standard error and Sigma for this specific run
        standard_error = std_BSM / sqrt(N_cycles);
        current_sigma = mean_BSM / standard_error; 
        
        sigma_values(counter) = current_sigma;
        
        % Print a progress update every 6 runs so you know it hasn't frozen
        if mod(counter, 6) == 0
            fprintf('  Completed %d / %d runs...\n', counter, total_iterations);
        end
        
        counter = counter + 1;
    end
end
fprintf('Monte Carlo Validation completed in %.2f seconds.\n', toc);

% --- NORMALITY ANALYSIS ---
mean_sigma = mean(sigma_values);
std_sigma = std(sigma_values);

fprintf('\n=== TRUE-SEED STATISTICAL CONFIDENCE ===\n');
fprintf('  Mean Extracted Sigma : %.2f \\sigma\n', mean_sigma);
fprintf('  Standard Deviation   : %.4f\n', std_sigma);

[h_ad, p_ad] = adtest(sigma_values);
if h_ad == 0
    fprintf('  Anderson-Darling Test: PASS (p = %.4f). Normal distribution confirmed.\n', p_ad);
else
    fprintf('  Anderson-Darling Test: FAIL (p = %.4f). Not normally distributed.\n', p_ad);
end

[h_lill, p_lill] = lillietest(sigma_values);
if h_lill == 0
    fprintf('  Lilliefors Test      : PASS (p = %.4f). Normal distribution confirmed.\n', p_lill);
else
    fprintf('  Lilliefors Test      : FAIL (p = %.4f). Not normally distributed.\n', p_lill);
end
fprintf('========================================\n');

%% ========================================================================
% PART 3: PUBLICATION-READY INDIVIDUAL VISUALIZATIONS
% =========================================================================
fprintf('\nGenerating Presentation Graphics (5 Separate Figures)...\n');

set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

color_Rb85 = [0 0.4470 0.7410];
color_Rb87 = [0.8500 0.3250 0.0980];
color_Signal = [0.4660 0.6740 0.1880];
color_Hist = [0.3010 0.7450 0.9330];

%%

% -------------------------------------------------------------------------
% FIGURE 1: Ab Initio Spatial Separation
% -------------------------------------------------------------------------
fig1 = figure('Name', 'Fig 1: Spatial Separation', 'Color', 'w', 'Position', [100, 100, 600, 450]);
hold on;
drop_center = -0.5 * g * T_drift^2;
zoom_idx = abs(x_vals - drop_center) < 60e-6; 

plot(x_vals(zoom_idx)*1e6, density_g_halfway(zoom_idx), 'Color', color_Rb85, 'LineWidth', 2.5);
plot(x_vals(zoom_idx)*1e6, density_e_halfway(zoom_idx), 'Color', color_Rb87, 'LineWidth', 2.5);
fill(x_vals(zoom_idx)*1e6, density_g_halfway(zoom_idx), color_Rb85, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
fill(x_vals(zoom_idx)*1e6, density_e_halfway(zoom_idx), color_Rb87, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

title('\textbf{Ab Initio Spatial Separation at $t=T$}', 'FontSize', 15);
xlabel('Vertical Position ($\mu$m)', 'FontSize', 14);
ylabel('Probability Density $|\psi(x)|^2$', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
legend('Lower Arm $|g\rangle$', 'Upper Arm $|e\rangle$ (Recoiled)', 'Location', 'northeast', 'FontSize', 11);

%%

% -------------------------------------------------------------------------
% FIGURE 2: Mid-Fringe Noise Rejection
% -------------------------------------------------------------------------
fig2 = figure('Name', 'Fig 2: Mid-Fringe Covariance', 'Color', 'w', 'Position', [150, 150, 500, 500]);
hold on;

% Visualization of mid-fringe covariance trace
scatter(demo_z85, demo_z87, 18, [0 0.1 0.5], 'filled', 'MarkerFaceAlpha', 0.9);

theta = linspace(0, 2*pi, 200);
x_ideal = -C_tau * sin(theta);
y_ideal = -C_tau * sin(theta + demo_dphi + pi/2); 
plot(x_ideal, y_ideal, 'Color', color_Rb87, 'LineWidth', 2);

axis square; grid on; box on;
title('\textbf{Mid-Fringe Phase Sweeping}', 'FontSize', 15);
xlabel('$^{85}$Rb Population Imbalance ($z_{85}$)', 'FontSize', 14);
ylabel('$^{87}$Rb Population Imbalance ($z_{87}$)', 'FontSize', 14);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlim([-1.2, 1.2]); ylim([-1.2, 1.2]);

% The legend string keeps the tildes to prevent PDF text clipping
legend('Active Sweep Measurements~~~~', 'Ideal Mid-Fringe Trace', 'Location', 'northeast', 'FontSize', 11);
%%

% -------------------------------------------------------------------------
% FIGURE 3: Lock-In Quadrature Histogram
% -------------------------------------------------------------------------
fig3 = figure('Name', 'Fig 3: Lock-In Extraction', 'Color', 'w', 'Position', [200, 200, 600, 450]);
hold on;

% Convert the data to microradians (urad)
data_urad = S_BSM_results * 1e6;
target_urad = target_signal * 1e6;

histogram(data_urad, 40, 'Normalization', 'pdf', 'FaceColor', color_Signal, 'EdgeColor', 'w', 'FaceAlpha', 0.85);
pd_phase = fitdist(data_urad, 'Normal');
x_fit_phase = linspace(min(data_urad), max(data_urad), 100);
plot(x_fit_phase, pdf(pd_phase, x_fit_phase), 'k-', 'LineWidth', 2.5);
xline(target_urad, 'r--', 'LineWidth', 3);

title('\textbf{Double-Difference Lock-In Isolation}', 'FontSize', 15);
xlabel('Isolated $S_{\rm BSM}$ Phase ($\mu$rad)', 'FontSize', 14);
ylabel('Probability Density', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% The legend string keeps the tildes to prevent PDF text clipping
legend('Post-Veto MC Data~~~~', 'Gaussian Fit', sprintf('Target: %.1f $\\mu$rad~~~~', target_urad), 'Location', 'northwest', 'FontSize', 11);
%%

% -------------------------------------------------------------------------
% FIGURE 4: Distribution of Sigma (Validation of Normality)
% -------------------------------------------------------------------------
fig4 = figure('Name', 'Fig 4: Sigma Distribution', 'Color', 'w', 'Position', [250, 250, 600, 450]);
hold on;

histogram(sigma_values, 10, 'Normalization', 'pdf', 'FaceColor', color_Hist, 'EdgeColor', 'w', 'FaceAlpha', 0.85);
pd_sigma = fitdist(sigma_values, 'Normal');
x_fit_sigma = linspace(min(sigma_values)-1, max(sigma_values)+1, 100);
plot(x_fit_sigma, pdf(pd_sigma, x_fit_sigma), 'k-', 'LineWidth', 2.5);
xline(mean_sigma, 'b--', 'LineWidth', 2.5);

title('\textbf{True-Seed Validation Ensembles}', 'FontSize', 15);
xlabel('Statistical Significance ($\sigma$)', 'FontSize', 14);
ylabel('Probability Density', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
legend(sprintf('%d Indep. MC Runs', length(sigma_values)), 'Gaussian Fit', sprintf('Mean: %.1f\\sigma', mean_sigma), 'Location', 'northeast', 'FontSize', 11);

%% 

% -------------------------------------------------------------------------
% FIGURE 5: Metrological Convergence (1/sqrt(N) scaling)
% -------------------------------------------------------------------------
fig5 = figure('Name', 'Fig 5: Integration Convergence', 'Color', 'w', 'Position', [300, 300, 600, 450]);
hold on;

% Calculate cumulative standard error
cycles_array = round(logspace(1, log10(N_cycles), 50));
std_err_array = zeros(size(cycles_array));
for i = 1:length(cycles_array)
    std_err_array(i) = std(S_BSM_results(1:cycles_array(i))) / sqrt(cycles_array(i)) * 1e9; % in nrad
end

% Plot actual data vs ideal 1/sqrt(N)
loglog(cycles_array, std_err_array, 'o-', 'Color', color_Signal, 'LineWidth', 2, 'MarkerFaceColor', color_Signal);
ideal_line = std_err_array(1) .* sqrt(cycles_array(1)) ./ sqrt(cycles_array);
loglog(cycles_array, ideal_line, 'k--', 'LineWidth', 2);

title('\textbf{Metrological Integration Convergence}', 'FontSize', 15);
xlabel('Number of Measurement Cycles ($N$)', 'FontSize', 14);
ylabel('Phase Standard Error (nrad)', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XScale', 'log', 'YScale', 'log');
legend('Simulated Phase Error', 'Ideal $1/\sqrt{N}$ SQL Limit', 'Location', 'southwest', 'FontSize', 11);

%%
% -------------------------------------------------------------------------
% FIGURE 6: Custom Q-Q Plot Normality Validation
% -------------------------------------------------------------------------
fig6 = figure('Name', 'Fig 6: Q-Q Plot', 'Color', 'w', 'Position', [400, 400, 500, 500]);
hold on;

% 1. Mathematical construction of the Q-Q data
N_samples = length(sigma_values);
sorted_data = sort(sigma_values);

% Calculate empirical cumulative probabilities using the standard (i - 0.5)/N formula
p_vals = ((1:N_samples)' - 0.5) / N_samples;

% Calculate theoretical quantiles for a Standard Normal distribution
theo_quantiles = norminv(p_vals, 0, 1);

% 2. Plot the data points
scatter(theo_quantiles, sorted_data, 50, [0 0.1 0.5], 'filled', ...
    'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.8);

% 3. Calculate and plot the ideal Gaussian reference line
% A perfect normal distribution follows the line: y = mu + sigma * x
sample_mu = mean(sigma_values);
sample_sig = std(sigma_values);

% Extend the dashed line slightly past the data points for visual clarity
x_line = linspace(min(theo_quantiles)-0.2, max(theo_quantiles)+0.2, 100);
y_line = sample_mu + sample_sig * x_line;

plot(x_line, y_line, 'k--', 'LineWidth', 2);

title('\textbf{Normality Validation (Q-Q Plot)}', 'FontSize', 15);
xlabel('Standard Normal Quantiles', 'FontSize', 14);
ylabel('Extracted Phase Quantiles ($\sigma$)', 'FontSize', 14);

grid on; box on; axis square;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% Safely apply the legend with the tilde trick
legend('120 MC Extractions~~~~', 'Ideal Gaussian Reference', 'Location', 'northwest', 'FontSize', 11);

%% ========================================================================
% PART 4: EXACT NUMERIC OUTPUTS FOR LATEX MANUSCRIPT
% =========================================================================
fprintf('\n=======================================================\n');
fprintf(' EXACT NUMERIC VALUES FOR LATEX MANUSCRIPT\n');
fprintf('=======================================================\n');
fprintf('TARGET SIGNAL      : %.1f nrad\n', target_signal * 1e9);
fprintf('EXTRACTED MEAN     : %.1f nrad\n', pd_phase.mu);
fprintf('STANDARD ERROR     : %.1f nrad\n', pd_phase.sigma/sqrt(N_cycles));
fprintf('CONFIDENCE (SIGMA) : %.2f σ\n', mean_sigma);
fprintf('P-VALUE (A-D TEST) : %.4f\n', p_ad);
fprintf('=======================================================\n\n');

%% ========================================================================
% LOCAL FUNCTIONS 
% =========================================================================
function [g_out, e_out] = apply_pulse(g_in, e_in, Area, k_vec, x)
    phase_fwd = exp(1i * k_vec * x);
    phase_bwd = exp(-1i * k_vec * x);
    C = cos(Area/2); S = sin(Area/2);
    g_out = C * g_in - 1i * S * phase_bwd .* e_in;
    e_out = C * e_in - 1i * S * phase_fwd .* g_in;
end

function [g_out, e_out] = apply_drift(g_in, e_in, U_T, U_V_half)
    g_out = U_V_half .* ifft(U_T .* fft( U_V_half .* g_in ));
    e_out = U_V_half .* ifft(U_T .* fft( U_V_half .* e_in ));
end

function delta_phi = extract_differential_phase(z1, z2, C_tau)
    % Unbiased Covariance Estimator for Mid-Fringe Operation
    % z1 = -A sin(phi)
    % z2 = -A sin(phi + delta_phi + pi/2) = -A cos(phi + delta_phi)
    % Cov(z1, z2) = - (A^2 / 2) * sin(delta_phi)
    
    C = cov(z1, z2);
    cov_z = C(1,2);
    
    sin_delta = - (2 * cov_z) / (C_tau^2);
    
    % Clamp strictly for mathematical safety at extreme noise bounds
    sin_delta = max(min(sin_delta, 1), -1); 
    delta_phi = asin(sin_delta);
end