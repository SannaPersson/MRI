
clear all; close all; clc;
%% Constants

gammabar = 42.58e6;
gamma = 2*pi*gammabar;
B0 = 3;

%% Hello world

iv = ImagingVolume(0, 0, 1e9, 1e9)

plot(iv)


%%

rf = RectPulse(5.9e-6, B0*gammabar, pi/2, 1e-3)


powspec(rf);


seemri(iv, B0, rf);

iv.toEquilibrium();

%% Experiments

% Original 
rf1 = RectPulse(5.9e-6, B0*gammabar, pi/2, 1e-3)

powspec(rf1);
%%

% Modify amplitude
rf2 = RectPulse(5.9e-7, B0*gammabar, pi/2, 1e-3)

powspec(rf2);

seemri(iv, B0, rf2);

iv.toEquilibrium();

% The amplitude changes how fast the signal changes

%% Relaxation
B0=3;
T1= 1.000;
T2 = 0.10;
iv=ImagingVolume(0,0,T1,T2);



rf = RectPulse(5.9e-7, B0*gammabar, pi/2, 0.2)

sampler=ADC(0.2, 0.2/100);

[S, ts] = seemri(iv, B0, rf, [], [], sampler);
size(S)
%plot(iv)
%%
size(S)
figure()
plot(ts, abs(S))
hold on;
plot(ts, angle(S))
legend(["Signal", "Phase"])



%% Relaxation: Run simulations
close all; clear all; clc;

% Constants
gammabar = 42.58e6;
gamma = 2*pi*gammabar;
B0 = 3;
B1 = 50e-6;
alpha = pi/2;
T1_values = [100e-3, 500e-3, 1000e-3]; % T1 values in seconds
T2_values = [10e-3, 50e-3, 100e-3]; % T2 values in seconds



% Initialize arrays to store simulation results
tp = alpha/ (gamma * B1) % pulse duration in seconds
rf = RectPulse(B1, B0*gammabar, 0, tp);
sampler = ADC(0.2, 0.2/100);

S_values = cell(length(T1_values), length(T2_values));
ts_values = cell(length(T1_values), length(T2_values));

% Run simulations
for i = 1:length(T1_values)
    for j = 1:length(T2_values)   
        rf = RectPulse(B1, B0*gammabar, pi/2, tp);
        sampler = ADC(0.2, 0.2/100);
        iv=ImagingVolume(0,0,T1_values(i),T2_values(j));

        [S, ts] = seemri(iv, B0, rf, [], [], sampler);
        iv.toEquilibrium();
        
        % Save results
        S_values{i, j} = S;
        ts_values{i, j} = ts;
    end
end


%% Plot results 
% Choose a T1, T2 combination for the phase plot
T1_index = 1; % change this to choose a different T1 value
T2_index = 1; % change this to choose a different T2 value

% Retrieve results for the chosen T1, T2 combination
S = S_values{T1_index, T2_index};
ts = ts_values{T1_index, T2_index};

% Create a figure for the signal and phase
figure;
subplot(2, 1, 1);
plot(ts, abs(S));
title(['T1 = ' num2str(T1_values(T1_index)) ' s, T2 = ' num2str(T2_values(T2_index)) ' s']);
xlabel('Time (s)');
ylabel('Magnitude');

subplot(2, 1, 2);
plot(ts, angle(S));
xlabel('Time (s)');
ylabel('Phase');

% Create a figure for the signals from all T1 and T2 combinations
figure;
for i = 1:length(T1_values)
    for j = 1:length(T2_values)
        % Retrieve results
        S = S_values{i, j};
        ts = ts_values{i, j};
        
        % Plot magnitude
        subplot(length(T1_values), length(T2_values), (i-1)*length(T2_values)+j);
        plot(ts, abs(S));
        title(['T1 = ' num2str(T1_values(i)) ' s, T2 = ' num2str(T2_values(j)) ' s']);
        xlabel('Time (s)');
        ylabel('Magnitude');
    end
end


%% Spin echo
clear all; close all; clc;
% Constants
T1 = 0.5; 
T2 = 20e-3; 
B0 = 3;
B1_90 = 50e-6;
gammabar = 42.58e6; 
gamma = 2*pi*gammabar; 
tau = 7e-2; % time spacing in seconds

alpha_1 = pi/2;
alpha_2 = pi;

tp = alpha_1/(gamma*B1_90)
B1_180 = alpha_2/(gamma*tp)

% Create the imaging volume with B0 inhomogeneity
iv = ImagingVolume(-4:1:4, -4:1:4, T1, T2, 1, 'dB0Gamma', 1e-6);

% Set the time spacing (tau) and calculate the pulse durations
B1 = [B1_90, B1_180, B1_180, B1_180, B1_180]
% Create the RF pulses for a CPMG sequence
theta_rf = [0, pi, pi, pi, pi]; % phases for the 90 and 180 degree pulses
%rf = RectPulse([B1,B1] gammabar*B0, theta_rf, tp, [0, tau]);
rf = RectPulse(B1, B0*gammabar, theta_rf, tp, [0, tau, 3*tau, 5*tau, 7*tau]);

sampler=ADC(0.05, 0.5e-3);
[S, ts] = seemri(iv, B0, rf, [], [], sampler);
iv.toEquilibrium();


%% Gradient echo
clear all; close all; clc;
% Constants
T1 = 0.5; 
T2 = 20e-3; 
B0 = 7;
B1 = 50e-6;
gammabar = 42.58e6; 
gamma = 2*pi*gammabar; 
tau = 2e-3; % time spacing in seconds
G = 1.2e-6
alpha_1 = pi/2;

tp = alpha_1/(gamma*B1)

% Create the imaging volume with B0 inhomogeneity
iv = ImagingVolume(-4:1:4, -4:1:4, T1, T2, 1, 'dB0Gamma', 1e-7);

rf = RectPulse(B1, B0*gammabar, 0, tp);

sampler=ADC(0.02, 0.02/100);
Gx = Gradient([tp tp+tau tp+3*tau], [-G 2*G 0]);
[S, ts] = seemri(iv, B0, rf, Gx, [], sampler);
iv.toEquilibrium();






%% Plot gradient pulse
[S_no, ts] = seemri(iv, B0, rf, [], [], sampler);
iv.toEquilibrium();
plot(ts, abs(S))
hold on;
plot(ts, abs(S_no))

legend(["Gradient echo", "Only RF"])



%% Slice selection
clear all; close all; clc;
% Constants
gammabar = 42.58e6; 
gamma = 2*pi*gammabar;
B1 = 2.93e-6;
tp = 2e-3;
alpha=gamma*B1*tp
B0 = 3;
iv1 = ImagingVolume(-1:0.25:5, -1:0.25:1, 0.8, 0.07, 1, 'PlotScale', 2);


rf = SincPulse(B1, gammabar*B0+8e3, 0, tp);
g = Gradient([0 tp 1.5*tp], 2*[70e-6 -70e-6 0]);
seemri(iv1, B0, rf, g, g, ADC(1.5*tp, tp/100));

