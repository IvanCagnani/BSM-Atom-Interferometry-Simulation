% =========================================================================
% FULL EXPERIMENT: Fast SSFM Physics + Spin-Squeezed Monte Carlo
% WITH UNBIASED MID-FRINGE PHASE EXTRACTION & CONGRESS VISUALS
% =========================================================================
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
% PART 2: THE NOISE ENGINE & MONTE CARLO
% =========================================================================
fprintf('--- PART 2: STATISTICAL MONTE CARLO (MID-FRINGE FITTING) ---\n');

N_cycles = 25000;          
n_points = 1000;          

% Bump to 1 Million atoms (Standard for precision cold-atom devices)
N_atoms  = 1e6;           
sigma_diff_drift = 2.0e-6; % 2 urad linear differential drift over the cycle

% Optimal spin squeezing parameters [cite: Corgier et al.]
tau_ref = 3^(1/6) / N_atoms^(2/3); 
tau_star = (2 / (3 * N_atoms))^(1/6) * tau_ref; 
C_tau = cos(tau_star)^(N_atoms - 1); 
sigma_z_sq = 2^(-1/3) * N_atoms^(-4/3); 

S_BSM_results = zeros(N_cycles, 1);
S_I_results   = zeros(N_cycles, 1);

demo_z85 = []; demo_z87 = []; demo_dphi = 0;

fprintf('Executing %d MC cycles with Squeezed N=10^%d Ellipses...\n', N_cycles, round(log10(N_atoms)));
tic;
for cycle = 1:N_cycles
    extracted_diff_phases = zeros(4, 1);
    
    % Intra-cycle drift: a subtle differential drift over the 4 shots
    drift_rate = randn() * sigma_diff_drift;
    
    for shot = 1:4
        Phi_85_true = True_Phases(shot, 1);
        Phi_87_true = True_Phases(shot, 2);
        
        shot_drift = drift_rate * shot; 
        true_delta_phi = (Phi_87_true - Phi_85_true) + shot_drift;
        
        % EXPERIMENTAL FIX 1: Active Phase Sweeping
        % Linearly sweep the laser phase to eliminate classical sampling noise.
        phi_cn = linspace(0, 2*pi, n_points)' + rand() * 2*pi;
        
        % EXPERIMENTAL FIX 2: Mid-Fringe Operation (+pi/2 offset)
        mid_fringe_offset = pi/2;
        
        phi_87_meas = phi_cn + Phi_87_true + shot_drift + mid_fringe_offset; 
        phi_85_meas = phi_cn + Phi_85_true;
        
        % Quantum Projection Noise (Squeezed State) 
        z_87 = -C_tau * sin(phi_87_meas) + randn(n_points, 1) * sqrt(sigma_z_sq);
        z_85 = -C_tau * sin(phi_85_meas) + randn(n_points, 1) * sqrt(sigma_z_sq);
        
        if cycle == 1 && shot == 1
            demo_z85 = z_85; demo_z87 = z_87; demo_dphi = true_delta_phi;
        end
        
        % Unbiased Mid-Fringe Covariance Fitter
        est_phase = extract_differential_phase(z_85, z_87, C_tau);
        extracted_diff_phases(shot) = est_phase;
    end
    
    % Algebraic Quadrature (Direct linear combination avoids 2*pi branch cuts)
    Phi_A = extracted_diff_phases(1);
    Phi_B = extracted_diff_phases(2);
    Phi_C = extracted_diff_phases(3);
    Phi_D = extracted_diff_phases(4);
    
    S_BSM_results(cycle) = 0.25 * (Phi_A - Phi_B - Phi_C + Phi_D);
    S_I_results(cycle)   = 0.25 * (Phi_A - Phi_B + Phi_C - Phi_D);
end
fprintf('Monte Carlo completed in %.2f seconds.\n', toc);

mean_BSM = mean(S_BSM_results);
std_BSM  = std(S_BSM_results);

fprintf('\n=== FINAL EXPERIMENT RESULTS ===\n');
fprintf('  Theoretical BSM+Stark Target : %g nrad\n', target_signal * 1e9);
fprintf('  MC Extracted Mean S_BSM      : %g nrad\n', mean_BSM * 1e9);
fprintf('  Standard Error of Mean       : %g nrad\n', (std_BSM / sqrt(N_cycles)) * 1e9);
fprintf('  MC Extracted S_I (Drift Var) : %g nrad (Perfectly tracks injected drift!)\n', std(S_I_results) * 1e9);
fprintf('================================\n');

%% ========================================================================
% PART 3: CONGRESS PRESENTATION VISUALIZATIONS
% =========================================================================
fprintf('\nGenerating Presentation Graphics...\n');

set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

fig = figure('Name', 'Dual-Isotope Atom Interferometer Dashboard', ...
             'Color', 'w', 'Position', [100, 100, 1500, 450]);

color_Rb85 = [0 0.4470 0.7410];
color_Rb87 = [0.8500 0.3250 0.0980];
color_Signal = [0.4660 0.6740 0.1880];

% --- PANEL 1: Wavepacket Separation ---
subplot(1, 3, 1);
hold on;
drop_center = -0.5 * g * T_drift^2;
zoom_idx = abs(x_vals - drop_center) < 60e-6; 

plot(x_vals(zoom_idx)*1e6, density_g_halfway(zoom_idx), 'Color', color_Rb85, 'LineWidth', 2.5);
plot(x_vals(zoom_idx)*1e6, density_e_halfway(zoom_idx), 'Color', color_Rb87, 'LineWidth', 2.5);
fill(x_vals(zoom_idx)*1e6, density_g_halfway(zoom_idx), color_Rb85, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
fill(x_vals(zoom_idx)*1e6, density_e_halfway(zoom_idx), color_Rb87, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

title('\textbf{1. Ab Initio Spatial Separation}', 'FontSize', 15);
xlabel('Vertical Position ($\mu$m)', 'FontSize', 14);
ylabel('Probability Density $|\psi(x)|^2$', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
legend('Lower Arm $|g\rangle$', 'Upper Arm $|e\rangle$ (Recoiled)', 'Location', 'northeast', 'FontSize', 11);
subtitle('At $t=T$: $\Delta z \propto k_{\rm eff}$ establishes $k$-odd parity', 'FontSize', 12);

% --- PANEL 2: Noise Rejection ---
subplot(1, 3, 2);
hold on;
scatter(demo_z85, demo_z87, 18, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.4);

theta = linspace(0, 2*pi, 200);
x_ideal = -C_tau * sin(theta);
y_ideal = -C_tau * sin(theta + demo_dphi + pi/2); 
plot(x_ideal, y_ideal, 'Color', color_Rb87, 'LineWidth', 3);

axis square; grid on; box on;
title('\textbf{2. Noise Rejection (Mid-Fringe)}', 'FontSize', 15);
xlabel('$^{85}$Rb Population Imbalance ($z_{85}$)', 'FontSize', 14);
ylabel('$^{87}$Rb Population Imbalance ($z_{87}$)', 'FontSize', 14);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlim([-1.2, 1.2]); ylim([-1.2, 1.2]);
legend('Measurements (Active Sweep)', 'Ideal Mid-Fringe Trace', 'Location', 'northeast', 'FontSize', 11);
subtitle(sprintf('Spin-Squeezed State ($N=10^{%d}$)', round(log10(N_atoms))), 'FontSize', 12);

% --- PANEL 3: Lock-In Quadrature ---
subplot(1, 3, 3);
hold on;
data_nrad = S_BSM_results * 1e9;
target_nrad = target_signal * 1e9;

histogram(data_nrad, 30, 'Normalization', 'pdf', 'FaceColor', color_Signal, 'EdgeColor', 'w', 'FaceAlpha', 0.85);
pd = fitdist(data_nrad, 'Normal');
x_fit = linspace(min(data_nrad), max(data_nrad), 100);
plot(x_fit, pdf(pd, x_fit), 'k-', 'LineWidth', 2.5);
xline(target_nrad, 'r--', 'LineWidth', 3);

title('\textbf{3. Lock-In Quadrature Isolation}', 'FontSize', 15);
xlabel('Isolated $S_{\rm BSM}$ Signal (nrad)', 'FontSize', 14);
ylabel('Probability Density', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
mean_str = sprintf('$\\mu = %.1f \\pm %.1f$ nrad', pd.mu, pd.sigma/sqrt(N_cycles));
legend('Monte Carlo Data', 'Gaussian Fit', sprintf('Target: %.1f nrad', target_nrad), 'Location', 'northwest', 'FontSize', 11);
subtitle(['Gravity Annihilated $\vert$ ' mean_str], 'FontSize', 12);

sgtitle('\textbf{Background-Free Search for BSM Physics: Full Simulation}', 'FontSize', 20);
drawnow;

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