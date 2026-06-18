%% Wind Turbine Aeroelasticity
%% TU Delft

clc
clear
close all

%% LOAD OPERATING CONDITIONS

load STATE

v0    = WindSpeeds;
omega = RtSpeeds*2*pi/60;
pitch = PitchAngles;

%% STRUCTURAL DATA

r = [1.50, 1.70, 2.70, 3.70, 4.70, 5.70, 6.70, 7.70, 8.70, 9.70, ...
     10.70, 11.70, 12.70, 13.70, 14.70, 15.70, 16.70, 17.70, ...
     19.70, 21.70, 23.70, 25.70, 27.70, 29.70, 31.70, 33.70, ...
     35.70, 37.70, 39.70, 41.70, 43.70, 45.70, 47.70, 49.70, ...
     51.70, 53.70, 55.70, 56.70, 57.70, 58.70, 59.20, 59.70, ...
     60.20, 60.70, 61.20, 61.70, 62.20, 62.70, 63.00];

m = [678.935, 678.935, 773.363, 740.550, 740.042, 592.496, ...
     450.275, 424.054, 400.638, 382.062, 399.655, 426.321, ...
     416.820, 406.186, 381.420, 352.822, 349.477, 346.538, ...
     339.333, 330.004, 321.990, 313.820, 294.734, 287.120, ...
     263.343, 253.207, 241.666, 220.638, 200.293, 179.404, ...
     165.094, 154.411, 138.935, 129.555, 107.264, 98.776, ...
     90.248, 83.001, 72.906, 68.772, 66.264, 59.340, ...
     55.914, 52.484, 49.114, 45.818, 41.669, 11.453, 10.319];

EI_flap = [18110.00E6, 18110.00E6, 19424.90E6, 17455.90E6, ...
           15287.40E6, 10782.40E6, 7229.72E6, 6309.54E6, ...
           5528.36E6, 4980.06E6, 4936.84E6, 4691.66E6, ...
           3949.46E6, 3386.52E6, 2933.74E6, 2568.96E6, ...
           2388.65E6, 2271.99E6, 2050.05E6, 1828.25E6, ...
           1588.71E6, 1361.93E6, 1102.38E6, 875.80E6, ...
           681.30E6, 534.72E6, 408.90E6, 314.54E6, ...
           238.63E6, 175.88E6, 126.01E6, 107.26E6, ...
           90.88E6, 76.31E6, 61.05E6, 49.48E6, ...
           39.36E6, 34.67E6, 30.41E6, 26.52E6, ...
           23.84E6, 19.63E6, 16.00E6, 12.83E6, ...
           10.08E6, 7.55E6, 4.60E6, 0.25E6, 0.17E6];

EI_edge = [18113.60E6, 18113.60E6, 19558.60E6, 19497.80E6, ...
           19788.80E6, 14858.50E6, 10220.60E6, 9144.70E6, ...
           8063.16E6, 6884.44E6, 7009.18E6, 7167.68E6, ...
           7271.66E6, 7081.70E6, 6244.53E6, 5048.96E6, ...
           4948.49E6, 4808.02E6, 4501.40E6, 4244.07E6, ...
           3995.28E6, 3750.76E6, 3447.14E6, 3139.07E6, ...
           2734.24E6, 2554.87E6, 2334.03E6, 1828.73E6, ...
           1584.10E6, 1323.36E6, 1183.68E6, 1020.16E6, ...
           797.81E6, 709.61E6, 518.19E6, 454.87E6, ...
           395.12E6, 353.72E6, 304.73E6, 281.42E6, ...
           261.71E6, 158.81E6, 137.88E6, 118.79E6, ...
           101.63E6, 85.07E6, 64.26E6, 6.61E6, 5.01E6];

%% AERO DATA
R_aero = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, ...
          28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, ...
          52.7500, 56.1667, 58.9000, 61.6333];
Twist_aero = [13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, ...
              7.795, 6.544, 5.361, 4.188, 3.125, 2.319, ...
              1.526, 0.863, 0.370, 0.106];

twist_r = interp1(R_aero, Twist_aero, r, 'linear', 'extrap');
%% STRUCTURAL PARAMETERS

R = max(r);
zeta = 0.00477465;



%% ==========================================================
%% MASS MATRIX
%% ==========================================================

Mf = 0;
Me = 0;

for j = 1:length(r)-1

    dr = r(j+1)-r(j);

    x = r(j)/R;

    phi_f = phif(x);

    phi_e = phie(x);

    Mf = Mf + m(j)*phi_f^2*dr;
    Me = Me + m(j)*phi_e^2*dr;

end

M = [Mf 0;
     0 Me];

%% ==========================================================
%% STIFFNESS MATRIX
%% ==========================================================

Kf = 0;
Ke = 0;

for j = 1:length(r)-1

    dr = r(j+1)-r(j);

    x = r(j)/R;

    d2phi_f = d2phif(x,R);

    d2phi_e = d2phie(x,R);

    Kf = Kf + EI_flap(j)*d2phi_f^2*dr;
    Ke = Ke + EI_edge(j)*d2phi_e^2*dr;

end

K = [Kf 0;
     0 Ke];

%% ==========================================================
%% DAMPING MATRIX
%% ==========================================================

Cff = 2*zeta*sqrt(K(1,1)*M(1,1));
Cee = 2*zeta*sqrt(K(2,2)*M(2,2));

C = [Cff 0;
     0 Cee];

%% ==========================================================
%% NATURAL FREQUENCIES
%% ==========================================================

omega_flap = sqrt(K(1,1)/M(1,1));
omega_edge = sqrt(K(2,2)/M(2,2));

f_flap = omega_flap/(2*pi);
f_edge = omega_edge/(2*pi);

fprintf('\nFlap frequency = %.3f Hz\n',f_flap);
fprintf('Edge frequency = %.3f Hz\n',f_edge);

%% STORAGE

Qf = zeros(length(v0),1);
Qe = zeros(length(v0),1);
P  = zeros(length(v0),1);

%% ==========================================================
%% LOOP OVER OPERATING POINTS
%% ==========================================================



for i = 1:length(v0)
    [Rx,FN,FT,P(i)] = BEM(v0(i),omega(i),pitch(i));
    
    % Interpolate BEM loads onto your structural grid points
    FN_interp = interp1(Rx,FN,r,'linear','extrap');
    FT_interp = interp1(Rx,FT,r,'linear','extrap');
    
    for j = 1:length(r)-1
        dr = r(j+1)-r(j);
        x = r(j)/R;
        
        % Calculate total orientation angle in radians (Twist + Pitch)
        theta = deg2rad(twist_r(j) + pitch(i));
        
        % STEP 2 COUPLING: Rotate global loads into local structural coordinates
        F_flap = FN_interp(j)*cos(theta) + FT_interp(j)*sin(theta);
        F_edge = -FN_interp(j)*sin(theta) + FT_interp(j)*cos(theta);
        
        phi_f = phif(x);
        phi_e = phie(x);
        
        % Compute modal generalized loads using properly coupled forces
        Qf(i) = Qf(i) + F_flap*phi_f*dr;
        Qe(i) = Qe(i) + F_edge*phi_e*dr;
    end
end

%% ==========================================================
%% POWER CURVE
%% ==========================================================

figure(2)

plot(v0,P,'bo-','LineWidth',1.5)

xlabel('Wind Speed (m/s)')
ylabel('Power (W)')
title('Power Curve')
grid on

%% ==========================================================
%% GENERALIZED FORCES
%% ==========================================================

figure(3)

plot(v0,Qf,'r-o','LineWidth',1.5)
hold on

plot(v0,Qe,'b-o','LineWidth',1.5)

xlabel('Wind Speed (m/s)')
ylabel('Generalized Force')
title('Modal Generalized Loads')

legend('Flapwise','Edgewise')
grid on

%% ==========================================================
%% STATIC MODAL DISPLACEMENTS
%% ==========================================================

u_flap = Qf./K(1,1);
u_edge = Qe./K(2,2);

figure(4)

plot(v0,u_flap,'r-o','LineWidth',1.5)
hold on

plot(v0,u_edge,'b-o','LineWidth',1.5)

xlabel('Wind Speed (m/s)')
ylabel('Static Deformation (m)')
title('Static deformation from 3m/s to 25 m/s')

legend('Flapwise','Edgewise')
grid on

%% ==========================================================
%% ROOT BENDING AND TIP DEFLECTION AT 12 m/s
%% ==========================================================

idx12 = 11; %Index of v=12m/s

% Tip Deflection

u_tip_flap = u_flap(idx12) * phif(1);
u_tip_edge = u_edge(idx12) * phie(1);

% Root Bending
[Rx, FN_12, FT_12, ~] = BEM(v0(idx12), omega(idx12), pitch(idx12));

M_root_flap = 0;
M_root_edge = 0;

twist_Rx = interp1(R_aero, Twist_aero, Rx, 'linear', 'extrap');

for j = 1:length(Rx)-1
    dr = Rx(j+1) - Rx(j);
    moment_arm = Rx(j); 
    
    % Get local section angle at this specific BEM node
    theta_Rx = deg2rad(twist_Rx(j) + pitch(idx12));
    
    % Rotate the forces for your moment calculations
    F_flap_12 = FN_12(j)*cos(theta_Rx) + FT_12(j)*sin(theta_Rx);
    F_edge_12 = -FN_12(j)*sin(theta_Rx) + FT_12(j)*cos(theta_Rx);
    
    % Force * distance * segment width
    M_root_flap = M_root_flap + F_flap_12 * moment_arm * dr;
    M_root_edge = M_root_edge + F_edge_12 * moment_arm * dr;
end

%% ==========================================================
%% TIME AND FREQUENCY DOMAIN PLOTS (OVERLAPPED) AT 12 M/S
%% ==========================================================

% 1. Setup ODE parameters for 12 m/s steady wind
idx12 = 11; % Index for v = 12 m/s
F_aero_12 = [Qf(idx12); Qe(idx12)]; % Constant modal aerodynamic forcing at 12 m/s

% 2. Solve Equation of Motion using ode45
% Equation: M*ddx + C*dx + K*x = F_aero_12
tspan = [0 60]; % Simulate for 60 seconds to capture transient and steady state
X0 = zeros(4,1); % Initial conditions [x_f; x_e; dx_f; dx_e] (starts from rest)

[t, X] = ode45(@(t,X) aeroelastic_ode(t, X, M, C, K, F_aero_12), tspan, X0);

% Extract generalized displacements from the solver
x_flap_time = X(:,1);
x_edge_time = X(:,2);

% 3. Calculate Tip Deflection over time
% Physical deflection: u(r,t) = x(t) * phi(r/R)
time_tip_flap = x_flap_time * phif(1);
time_tip_edge = x_edge_time * phie(1);

% 4. Calculate Root Bending Moment over time
% Structural beam theory: M(r,t) = EI(r) * d^2u/dr^2
% Evaluated at root (x=0)
time_mom_flap = EI_flap(1) * x_flap_time * d2phif(0, R);
time_mom_edge = EI_edge(1) * x_edge_time * d2phie(0, R);

%% --- FIGURE 1: TIME DOMAIN PERFORMANCE (OVERLAPPED) ---
figure('Name', 'Time Domain Responses at 12 m/s', 'NumberTitle', 'off');

% Top Subplot: Tip Deflection Comparison
subplot(2,1,1);
plot(t, time_tip_flap, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, time_tip_edge, 'b-', 'LineWidth', 1.5); 
xlabel('Time (s)'); ylabel('Deflection (m)');
title('Tip Deflection at 12 m/s (Time Domain)');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

% Bottom Subplot: Root Bending Moment Comparison
subplot(2,1,2);
plot(t, time_mom_flap/1e6, 'r-', 'LineWidth', 1.5); % Convert to MNm
hold on;
plot(t, time_mom_edge/1e6, 'b-', 'LineWidth', 1.5); % Convert to MNm
xlabel('Time (s)'); ylabel('Moment (MNm)');
title('Root Bending Moment at 12 m/s (Time Domain)');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;


%% --- FIGURE 2: FREQUENCY DOMAIN PERFORMANCE (OVERLAPPED) ---
% Note: ode45 uses variable time steps. We must interpolate to a uniform 
% time grid before applying the Fast Fourier Transform (FFT).

t_unif = linspace(min(t), max(t), 2000); % 2000 points uniform grid
Fs = 1/mean(diff(t_unif)); % Sampling frequency
L = length(t_unif);
freq = Fs * (0:(L/2))/L; % Frequency axis

% Interpolate time domain signals
u_flap_unif = interp1(t, time_tip_flap, t_unif);
u_edge_unif = interp1(t, time_tip_edge, t_unif);
m_flap_unif = interp1(t, time_mom_flap/1e6, t_unif);
m_edge_unif = interp1(t, time_mom_edge/1e6, t_unif);

% Calculate FFTs (Single-Sided Amplitude Spectrum)
fft_u_flap = compute_fft(u_flap_unif, L);
fft_u_edge = compute_fft(u_edge_unif, L);
fft_m_flap = compute_fft(m_flap_unif, L);
fft_m_edge = compute_fft(m_edge_unif, L);

figure('Name', 'Frequency Domain Responses at 12 m/s', 'NumberTitle', 'off');

% Top Subplot: Deflection Frequency Comparison
subplot(2,1,1);
plot(freq, fft_u_flap, 'r-', 'LineWidth', 1.5);
hold on;
plot(freq, fft_u_edge, 'b-', 'LineWidth', 1.5); 
xlim([0, 3]); xlabel('Frequency (Hz)'); ylabel('Magnitude (m)');
title('Tip Deflection FFT Spectrum');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

% Bottom Subplot: Root Moment Frequency Comparison
subplot(2,1,2);
plot(freq, fft_m_flap, 'r-', 'LineWidth', 1.5);
hold on;
plot(freq, fft_m_edge, 'b-', 'LineWidth', 1.5);
xlim([0, 3]); xlabel('Frequency (Hz)'); ylabel('Magnitude (MNm)');
title('Root Bending Moment FFT Spectrum');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

%% ==========================================================
%% TASK 4: PERIODIC WIND CONDITIONS
%% ==========================================================

% 1. Parameters for periodic wind around 15 m/s
omega_fixed = 1.267; % rad/s (from Table 2 for 15 m/s)
pitch_fixed = 10.45; % degrees (from Table 1 for 15 m/s)

% 2. Create a Lookup Table for Aerodynamic Forces (14 m/s to 16 m/s)
% This avoids running the heavy BEM code at every ODE time step.
V_lookup = linspace(14, 16, 15);
Qf_lookup = zeros(size(V_lookup));
Qe_lookup = zeros(size(V_lookup));

fprintf('\nPre-calculating BEM lookup table for periodic wind...\n');
for k = 1:length(V_lookup)
    [Rx_per, FN_per, FT_per, ~] = BEM(V_lookup(k), omega_fixed, pitch_fixed);
    
    FN_interp = interp1(Rx_per, FN_per, r, 'linear', 'extrap');
    FT_interp = interp1(Rx_per, FT_per, r, 'linear', 'extrap');
    
    Qf_k = 0; Qe_k = 0;
    for j = 1:length(r)-1
        dr = r(j+1) - r(j);
        x = r(j)/R;
        theta = deg2rad(twist_r(j) + pitch_fixed);
        
        F_flap = FN_interp(j)*cos(theta) + FT_interp(j)*sin(theta);
        F_edge = -FN_interp(j)*sin(theta) + FT_interp(j)*cos(theta);
        
        Qf_k = Qf_k + F_flap * phif(x) * dr;
        Qe_k = Qe_k + F_edge * phie(x) * dr;
    end
    Qf_lookup(k) = Qf_k;
    Qe_lookup(k) = Qe_k;
end
fprintf('Lookup table completed.\n');

% 3. Solve ODE for Periodic Wind
tspan_per = [0 60]; % Simulate for 60 seconds
X0_per = zeros(4,1); % Start from rest [x_f; x_e; dx_f; dx_e]

[t_per, X_per] = ode45(@(t,X) aeroelastic_ode_periodic(t, X, M, C, K, V_lookup, Qf_lookup, Qe_lookup), tspan_per, X0_per);

x_flap_per = X_per(:,1);
x_edge_per = X_per(:,2);

% 4. Calculate Time-Domain Responses (Tip Deflection & Root Moment)
tip_flap_per = x_flap_per * phif(1);
tip_edge_per = x_edge_per * phie(1);

mom_flap_per = EI_flap(1) * x_flap_per * d2phif(0, R);
mom_edge_per = EI_edge(1) * x_edge_per * d2phie(0, R);

%% --- FIGURE 3: TIME DOMAIN PERFORMANCE (PERIODIC) ---
figure('Name', 'Time Domain Responses (Periodic Wind)', 'NumberTitle', 'off');

subplot(2,1,1);
plot(t_per, tip_flap_per, 'r-', 'LineWidth', 1.5); hold on;
plot(t_per, tip_edge_per, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Deflection (m)');
title('Tip Deflection under Periodic Wind');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t_per, mom_flap_per/1e6, 'r-', 'LineWidth', 1.5); hold on;
plot(t_per, mom_edge_per/1e6, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Moment (MNm)');
title('Root Bending Moment under Periodic Wind');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

%% --- FIGURE 4: FREQUENCY DOMAIN PERFORMANCE (PERIODIC) ---
% Interpolate to uniform time grid for FFT
t_unif_per = linspace(min(t_per), max(t_per), 3000); 
Fs_per = 1/mean(diff(t_unif_per)); 
L_per = length(t_unif_per);
freq_per = Fs_per * (0:(L_per/2))/L_per; 

u_flap_unif_per = interp1(t_per, tip_flap_per, t_unif_per);
u_edge_unif_per = interp1(t_per, tip_edge_per, t_unif_per);
m_flap_unif_per = interp1(t_per, mom_flap_per/1e6, t_unif_per);
m_edge_unif_per = interp1(t_per, mom_edge_per/1e6, t_unif_per);

fft_u_flap_per = compute_fft(u_flap_unif_per, L_per);
fft_u_edge_per = compute_fft(u_edge_unif_per, L_per);
fft_m_flap_per = compute_fft(m_flap_unif_per, L_per);
fft_m_edge_per = compute_fft(m_edge_unif_per, L_per);

figure('Name', 'Frequency Domain Responses (Periodic Wind)', 'NumberTitle', 'off');

subplot(2,1,1);
plot(freq_per, fft_u_flap_per, 'r-', 'LineWidth', 1.5); hold on;
plot(freq_per, fft_u_edge_per, 'b-', 'LineWidth', 1.5); 
xlim([0, 1.5]); xlabel('Frequency (Hz)'); ylabel('Magnitude (m)');
title('Tip Deflection FFT Spectrum (Periodic Wind)');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(freq_per, fft_m_flap_per, 'r-', 'LineWidth', 1.5); hold on;
plot(freq_per, fft_m_edge_per, 'b-', 'LineWidth', 1.5);
xlim([0, 1.5]); xlabel('Frequency (Hz)'); ylabel('Magnitude (MNm)');
title('Root Bending Moment FFT Spectrum (Periodic Wind)');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

%% ==========================================================
%% TASK 5: IMPACT OF BLADE DEFLECTIONS (AERODYNAMIC DAMPING)
%% ==========================================================

fprintf('Simulating Task 5 (with Aerodynamic Damping)...\n');

% Solve ODE activating the blade velocity term
[t_t5, X_t5] = ode45(@(t,X) aeroelastic_ode_task5(t, X, M, C, K, V_lookup, Qf_lookup, Qe_lookup), tspan_per, X0_per);

x_flap_t5 = X_t5(:,1);
x_edge_t5 = X_t5(:,2);

% Calculate Tip Deflection
tip_flap_t5 = x_flap_t5 * phif(1);
tip_edge_t5 = x_edge_t5 * phie(1);

% Calculate Root Bending Moment
mom_flap_t5 = EI_flap(1) * x_flap_t5 * d2phif(0, R);
mom_edge_t5 = EI_edge(1) * x_edge_t5 * d2phie(0, R);

%% --- FIGURE 5: COMPARISON (DEACTIVATED VS ACTIVATED) - SEPARATED MODES ---
% Create a larger figure to accommodate a 2x2 layout
figure('Name', 'Task 5: Aerodynamic Damping Comparison', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 700]);

% 1. Top-Left: Flapwise Tip Deflection
subplot(2,2,1);
plot(t_per, tip_flap_per, 'r--', 'LineWidth', 1.5); hold on; % Task 4 (Deactivated)
plot(t_t5, tip_flap_t5, 'b-', 'LineWidth', 1.5);             % Task 5 (Activated)
xlabel('Time (s)'); ylabel('Deflection (m)');
title('Flapwise Tip Deflection');
legend('Deactivated (Task 4)', 'Activated (Task 5)', 'Location', 'best');
grid on;

% 2. Top-Right: Edgewise Tip Deflection
subplot(2,2,2);
plot(t_per, tip_edge_per, 'r--', 'LineWidth', 1.5); hold on; 
plot(t_t5, tip_edge_t5, 'b-', 'LineWidth', 1.5);             
xlabel('Time (s)'); ylabel('Deflection (m)');
title('Edgewise Tip Deflection');
legend('Deactivated (Task 4)', 'Activated (Task 5)', 'Location', 'best');
grid on;

% 3. Bottom-Left: Flapwise Root Bending Moment
subplot(2,2,3);
plot(t_per, mom_flap_per/1e6, 'r--', 'LineWidth', 1.5); hold on; 
plot(t_t5, mom_flap_t5/1e6, 'b-', 'LineWidth', 1.5);             
xlabel('Time (s)'); ylabel('Moment (MNm)');
title('Flapwise Root Bending Moment');
grid on;

% 4. Bottom-Right: Edgewise Root Bending Moment
subplot(2,2,4);
plot(t_per, mom_edge_per/1e6, 'r--', 'LineWidth', 1.5); hold on; 
plot(t_t5, mom_edge_t5/1e6, 'b-', 'LineWidth', 1.5);             
xlabel('Time (s)'); ylabel('Moment (MNm)');
title('Edgewise Root Bending Moment');
grid on;

%% ==========================================================
%% TASK 6: CENTRIFUGAL AND GRAVITY STIFFENING
%% ==========================================================

fprintf('Calculating Geometric Stiffness for Task 6...\n');

g_gravity = 9.81; % Gravity (m/s^2)

% 1. Calculate Axial Tension distributions T_c(r) and T_g_amp(r)
% We integrate the forces from the tip to the root (backwards loop)
Tc = zeros(size(r));
Tg_amp = zeros(size(r));

for i = length(r)-1:-1:1
    dr = r(i+1) - r(i);
    r_mid = (r(i+1) + r(i))/2; % Center of the segment
    
    % Force element: mass/length * dr * acceleration
    dF_c = m(i) * omega_fixed^2 * r_mid * dr; % Centrifugal
    dF_g = m(i) * g_gravity * dr;             % Gravity amplitude
    
    % Tension accumulates towards the root
    Tc(i) = Tc(i+1) + dF_c;
    Tg_amp(i) = Tg_amp(i+1) + dF_g;
end

% 2. Calculate Modal Geometric Stiffness base terms
K_geo_c_f = 0; K_geo_c_e = 0;
K_geo_g_f = 0; K_geo_g_e = 0;

for j = 1:length(r)-1
    dr = r(j+1) - r(j);
    x = r(j)/R;
    
    df_f = dphif(x, R);
    df_e = dphie(x, R);
    
    % Centrifugal stiffening integrals
    K_geo_c_f = K_geo_c_f + Tc(j) * df_f^2 * dr;
    K_geo_c_e = K_geo_c_e + Tc(j) * df_e^2 * dr;
    
    % Gravity stiffening integrals
    K_geo_g_f = K_geo_g_f + Tg_amp(j) * df_f^2 * dr;
    K_geo_g_e = K_geo_g_e + Tg_amp(j) * df_e^2 * dr;
end

% 3. Solve ODE with Geometric Stiffness Active
fprintf('Simulating Task 6 (with Centrifugal & Gravity effects)...\n');

[t_t6, X_t6] = ode45(@(t,X) aeroelastic_ode_task6(t, X, M, C, K, ...
    K_geo_c_f, K_geo_c_e, K_geo_g_f, K_geo_g_e, omega_fixed, ...
    V_lookup, Qf_lookup, Qe_lookup), tspan_per, X0_per);
    
x_flap_t6 = X_t6(:,1);
tip_flap_t6 = x_flap_t6 * phif(1);

%% --- FIGURE 6: TASK 6 COMPARISON ---
figure('Name', 'Task 6: Geometric Stiffening', 'NumberTitle', 'off');

% Top Subplot: Time Domain
subplot(2,1,1);
plot(t_t5, tip_flap_t5, 'r--', 'LineWidth', 1.5); hold on; % Task 5 (No Geo Stiffening)
plot(t_t6, tip_flap_t6, 'b-', 'LineWidth', 1.5);           % Task 6 (With Geo Stiffening)
xlabel('Time (s)'); ylabel('Flapwise Deflection (m)');
title('Flapwise Tip Deflection: Geometric Stiffening (Time Domain)');
legend('Without Geo Stiffening (Task 5)', 'With Geo Stiffening (Task 6)', 'Location', 'best');
grid on;

% Bottom Subplot: Frequency Domain
% Interpolate for FFT
u_flap_unif_t6 = interp1(t_t6, tip_flap_t6, t_unif_per);
fft_u_flap_t6 = compute_fft(u_flap_unif_t6, L_per);

subplot(2,1,2);
plot(freq_per, fft_u_flap_per, 'r--', 'LineWidth', 1.5); hold on;
plot(freq_per, fft_u_flap_t6, 'b-', 'LineWidth', 1.5);
xlim([0, 2]); xlabel('Frequency (Hz)'); ylabel('Magnitude (m)');
title('Flapwise Tip Deflection FFT: Geometric Stiffening (Freq Domain)');
legend('Without Geo Stiffening', 'With Geo Stiffening', 'Location', 'best');
grid on;

%% ==========================================================
%% TASK 7: UNSTEADY LOADS (DYNAMIC INFLOW) - OPTION 1
%% ==========================================================

fprintf('Simulating Task 7 (Dynamic Inflow on Rotor Thrust)...\n');

v0_t7 = 15;        % Wind speed (m/s)
omega_t7 = 1.267;  % Rotor speed (rad/s)
f_pitch = [0.05, 0.2, 0.5]; % Frequencies (Hz)

% Time vector (simulate 60 seconds with a small time step for ODE accuracy)
dt = 0.05;
t_t7 = 0:dt:60; 

% Pre-allocate storage for Thrust
Thrust_QS = zeros(length(f_pitch), length(t_t7)); % Quasi-Steady (Without model)
Thrust_DI = zeros(length(f_pitch), length(t_t7)); % Dynamic Inflow (With model)

figure('Name', 'Task 7: Unsteady Rotor Loads', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 800]);

for k = 1:length(f_pitch)
    f = f_pitch(k);
    
    % Harmonic collective pitch
    pitch_t7 = 10.45 + 5 * sin(2 * pi * f * t_t7);
    
    % Initialize induction factors for the Dynamic Inflow time-marching
    a_dyn_prev = zeros(17, 1); % 17 aerodynamic nodes
    
    for i = 1:length(t_t7)
        current_pitch = pitch_t7(i);
        
        % 1. Quasi-Steady Thrust (Standard BEM - Instantaneous)
        % Note: We capture Rx_aero to integrate over the correct aerodynamic grid
        [Rx_aero, FN_qs, ~, ~] = BEM(v0_t7, omega_t7, current_pitch);
        
        % Integrate normal force along the 3 blades to get Rotor Thrust
        T_qs = 0;
        for j = 1:length(Rx_aero)-1
            T_qs = T_qs + FN_qs(j) * (Rx_aero(j+1) - Rx_aero(j));
        end
        Thrust_QS(k, i) = T_qs * 3; % 3 blades
        
        % 2. Dynamic Inflow Thrust (Time-filtered Induction)
        [~, FN_di, ~, ~, a_dyn_new] = BEM_dynamic(v0_t7, omega_t7, current_pitch, dt, a_dyn_prev);
        
        T_di = 0;
        for j = 1:length(Rx_aero)-1
            T_di = T_di + FN_di(j) * (Rx_aero(j+1) - Rx_aero(j));
        end
        Thrust_DI(k, i) = T_di * 3; 
        
        % Update induction state for the next time step
        a_dyn_prev = a_dyn_new;
    end
    
    % Plotting results for this frequency
    subplot(3, 1, k);
    plot(t_t7, Thrust_QS(k, :)/1000, 'r--', 'LineWidth', 1.5); hold on;
    plot(t_t7, Thrust_DI(k, :)/1000, 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Rotor Thrust (kN)');
    title(sprintf('Rotor Thrust Response at Pitch Frequency f = %.2f Hz', f));
    legend('Quasi-Steady (No Dynamic Inflow)', 'Unsteady (With Dynamic Inflow)', 'Location', 'best');
    grid on;
end

%% ==========================================================
%% HELPER FUNCTIONS
%% ==========================================================

function [phif] = phif(x)
         phif = ...
              0.0622*x^2 ...
            + 1.7254*x^3 ...
            - 3.2452*x^4 ...
            + 4.7131*x^5 ...
            - 2.2555*x^6;
end

function [phie] = phie(x)
        phie = ...
              0.3627*x^2 ...
            + 2.5337*x^3 ...
            - 3.5772*x^4 ...
            + 2.3760*x^5 ...
            - 0.6952*x^6;
end

function [d2phif] = d2phif(x, R)
    d2phif = ...
       (2*0.0622 ...
      +6*1.7254*x ...
      -12*3.2452*x^2 ...
      +20*4.7131*x^3 ...
      -30*2.2555*x^4)/R^2;
end

function [d2phie] = d2phie(x, R)
    d2phie = ...
       (2*0.3627 ...
      +6*2.5337*x ...
      -12*3.5772*x^2 ...
      +20*2.3760*x^3 ...
      -30*0.6952*x^4)/R^2;
end

function dXdt = aeroelastic_ode(t, X, M, C, K, F)
    % State-space formulation for 2-DOF Structural System
    % X = [x_flap; x_edge; dx_flap; dx_edge]
    
    x  = X(1:2); % Displacements
    dx = X(3:4); % Velocities
    
    % M*ddx + C*dx + K*x = F  =>  ddx = M \ (F - C*dx - K*x)
    ddx = M \ (F - C*dx - K*x);
    
    % Return state derivatives [velocities; accelerations]
    dXdt = [dx; ddx];
end

function P1 = compute_fft(signal, L)
    % Computes the single-sided amplitude spectrum of a signal
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
end

function dXdt = aeroelastic_ode_periodic(t, X, M, C, K, V_lookup, Qf_lookup, Qe_lookup)
    % 1. Calculate periodic wind speed at time t
    Vt = 15 + 0.5*cos(1.267*t) + 0.085*cos(2.534*t) + 0.015*cos(3.801*t);
    
    % 2. Interpolate generalized forces for current wind speed
    F_flap = interp1(V_lookup, Qf_lookup, Vt, 'linear', 'extrap');
    F_edge = interp1(V_lookup, Qe_lookup, Vt, 'linear', 'extrap');
    F = [F_flap; F_edge];
    
    % 3. State-space matrices
    x  = X(1:2);
    dx = X(3:4);
    
    % 4. Solve for accelerations
    ddx = M \ (F - C*dx - K*x);
    dXdt = [dx; ddx];
end

function dXdt = aeroelastic_ode_task5(t, X, M, C, K, V_lookup, Qf_lookup, Qe_lookup)
    % 1. Calculate periodic wind speed at time t
    Vt = 15 + 0.5*cos(1.267*t) + 0.085*cos(2.534*t) + 0.015*cos(3.801*t);
    
    x  = X(1:2);
    dx = X(3:4);
    
    % ===================================================================
    % TASK 5: ACTIVATE BLADE VELOCITY EFFECT (AERODYNAMIC DAMPING)
    % Eq 1: V_outplan = V0 - V_blade_out
    % ===================================================================
    
    % Calculate the mode shape value at a representative aerodynamic 
    % station (75% span) to estimate the effective blade velocity.
    x_75 = 0.75;
    phi_f_75 = 0.0622*x_75^2 + 1.7254*x_75^3 - 3.2452*x_75^4 + 4.7131*x_75^5 - 2.2555*x_75^6;
    
    % Out-of-plane blade velocity (Flapwise derivative)
    v_blade_out = dx(1) * phi_f_75; 
    
    % Effective relative wind speed seen by the blade
    V_eff = Vt - v_blade_out; 
    
    % Interpolate generalized forces for the EFFECTIVE wind speed
    F_flap = interp1(V_lookup, Qf_lookup, V_eff, 'linear', 'extrap');
    F_edge = interp1(V_lookup, Qe_lookup, V_eff, 'linear', 'extrap');
    
    F = [F_flap; F_edge];
    
    % Solve for accelerations
    ddx = M \ (F - C*dx - K*x);
    dXdt = [dx; ddx];
end

function [dphif] = dphif(x, R)
    % First spatial derivative of flapwise mode shape (d(phi_f)/dr)
    % Note: since x = r/R, d(phi)/dr = (1/R) * d(phi)/dx
    dphif = (2*0.0622*x + 3*1.7254*x^2 - 4*3.2452*x^3 + 5*4.7131*x^4 - 6*2.2555*x^5) / R;
end

function [dphie] = dphie(x, R)
    % First spatial derivative of edgewise mode shape
    dphie = (2*0.3627*x + 3*2.5337*x^2 - 4*3.5772*x^3 + 5*2.3760*x^4 - 6*0.6952*x^5) / R;
end

function dXdt = aeroelastic_ode_task6(t, X, M, C, K, Kc_f, Kc_e, Kg_f, Kg_e, omega, V_lookup, Qf_lookup, Qe_lookup)
    % Periodic wind
    Vt = 15 + 0.5*cos(1.267*t) + 0.085*cos(2.534*t) + 0.015*cos(3.801*t);
    
    x  = X(1:2);
    dx = X(3:4);
    
    % Aerodynamic Damping (from Task 5)
    x_75 = 0.75;
    phi_f_75 = 0.0622*x_75^2 + 1.7254*x_75^3 - 3.2452*x_75^4 + 4.7131*x_75^5 - 2.2555*x_75^6;
    v_blade_out = dx(1) * phi_f_75; 
    V_eff = Vt - v_blade_out; 
    
    F_flap = interp1(V_lookup, Qf_lookup, V_eff, 'linear', 'extrap');
    F_edge = interp1(V_lookup, Qe_lookup, V_eff, 'linear', 'extrap');
    F = [F_flap; F_edge];
    
    % ===================================================================
    % TASK 6: GEOMETRIC STIFFNESS MATRIX K_geo(t)
    % ===================================================================
    % The gravity term fluctuates at 1P (omega * t). It alternates between 
    % tension (increasing stiffness) and compression (decreasing stiffness).
    K_geo_flap = Kc_f - Kg_f * cos(omega * t);
    K_geo_edge = Kc_e - Kg_e * cos(omega * t);
    
    % Add geometric stiffness to the base structural stiffness
    K_total = K + diag([K_geo_flap, K_geo_edge]);
    
    % Solve for accelerations using the total stiffness
    ddx = M \ (F - C*dx - K_total*x);
    dXdt = [dx; ddx];
end