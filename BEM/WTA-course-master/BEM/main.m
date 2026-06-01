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

%% FUNCTION DEFINITION

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

% 1. Create a dummy time vector (0 to 10 seconds)
t = 0:0.1:10; 

% 2. Expand single static numbers into constant arrays over time
time_tip_flap = u_tip_flap * ones(size(t));
time_tip_edge = u_tip_edge * ones(size(t));
time_mom_flap  = M_root_flap * ones(size(t));
time_mom_edge  = M_root_edge * ones(size(t));

%% --- FIGURE 1: TIME DOMAIN PERFORMANCE (OVERLAPPED) ---
figure('Name', 'Time Domain Responses at 12 m/s', 'NumberTitle', 'off');

% Top Subplot: Tip Deflection Comparison
subplot(2,1,1);
plot(t, time_tip_flap, 'r-', 'LineWidth', 2);
hold on;
plot(t, time_tip_edge, 'b--', 'LineWidth', 2); % Dashed line for edgewise
xlabel('Time (s)'); ylabel('Deflection (m)');
title('Tip Deflection Comparison');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

% Bottom Subplot: Root Bending Moment Comparison
subplot(2,1,2);
plot(t, time_mom_flap/1e6, 'r-', 'LineWidth', 2); % Converted to MNm
hold on;
plot(t, time_mom_edge/1e6, 'b--', 'LineWidth', 2); % Converted to MNm
xlabel('Time (s)'); ylabel('Moment (MNm)');
title('Root Bending Moment Comparison');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;


%% --- FIGURE 2: FREQUENCY DOMAIN PERFORMANCE (OVERLAPPED) ---
figure('Name', 'Frequency Domain Responses at 12 m/s', 'NumberTitle', 'off');

% Setup a basic frequency axis
freq = 0:0.1:5; 

fft_tip_flap = zeros(size(freq)); fft_tip_flap(freq == 0) = u_tip_flap;
fft_tip_edge = zeros(size(freq)); fft_tip_edge(freq == 0) = u_tip_edge;
fft_mom_flap = zeros(size(freq)); fft_mom_flap(freq == 0) = M_root_flap/1e6;
fft_mom_edge = zeros(size(freq)); fft_mom_edge(freq == 0) = M_root_edge/1e6;

% Top Subplot: Deflection Frequency Comparison
subplot(2,1,1);
plot(freq, fft_tip_flap, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;
plot(freq, fft_tip_edge, 'b--x', 'LineWidth', 2, 'MarkerSize', 8); 
xlim([-0.5, 3]); xlabel('Frequency (Hz)'); ylabel('Magnitude (m)');
title('Tip Deflection FFT Magnitude');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

% Bottom Subplot: Root Moment Frequency Comparison
subplot(2,1,2);
plot(freq, fft_mom_flap, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;
plot(freq, fft_mom_edge, 'b--x', 'LineWidth', 2, 'MarkerSize', 8);
xlim([-0.5, 3]); xlabel('Frequency (Hz)'); ylabel('Magnitude (MNm)');
title('Root Bending Moment FFT Magnitude');
legend('Flapwise', 'Edgewise', 'Location', 'best');
grid on;

