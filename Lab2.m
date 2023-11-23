clear all; close all; clc;

gammabar = 42.58e6;
gamma = 2*pi*gammabar;
B0 = 1.5;



f_rf = gamma*B0;
tp = 1e-3;
alpha=pi/2;
B1 = alpha/(gamma*tp); 


rf = RectPulse(B1, f_rf, 0, tp);

FOV = 8e-3;
pixel_size = 1e-3; 

iv = disc(6, 0.5); % rad, res


% Given parameters
TE = 10e-3; 
TR = 2; 

kmaxx = 1/pixel_size; 
kmaxy = 1/pixel_size; 

dkx = 1/FOV;
dky = 1/FOV;

kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid

tau = 4e-3;
%Q4: We need  16 excitation repetitions, to measure each line in y directions
N = length(kya)
% Calculate the gradients as the delta_k/(gamma*tau)
Gpexs = kxa/(gamma*tau); 

gx = Gradient([tp tp+tau], {Gpexs 0});


% Frequency encoding gradients
Gfey1 = -2*(kmaxy/(gamma*tau));
Gfey2 = kmaxy/(gamma*tau);
% same area of both freq encoding pulses
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

% Need to sample twice times the frequency(Nyqvist)
max_freq = 2/pixel_size

dt = 1/max_freq

sampler = ADC(TE-tau, TE+tau, dt);

[S, ts] = seemri(iv, B0, rf, gx, gy, sampler, TR, length(Gpexs), 'PlotKSpace', true);

figure()

[imout, x, y] = mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true);

%% True image
figure()
side = sqrt(length(iv.Mz0))
side
imagesc(reshape(iv.Mz0, side, side))


%% Sketch of kspace
figure; 
hold on; 
kya_new = zeros(size(kya));
order = [8, 9, 7, 10, 6, 11, 5, 12, 4, 13, 3, 14, 2, 15, 1, 16];

for i = 1:length(order)
    plot(kxa, ones(size(kxa))*kya(order(i)), 'b-');
    
    quiver(kxa(1), kya(order(i)), kxa(2)-kxa(1), 0, 'r', 'MaxHeadSize', 0.05, 'AutoScale', 'off');
    
    text(kxa(1), kya(order(i)), num2str(i), 'VerticalAlignment', 'bottom');
end
xlabel('kx');
ylabel('ky');
title('k-space trajectories'); 

hold off; 


%% Design of sequence parameters
gammabar = 42.58e6;
gamma = 2*pi*gammabar;
B0 = 1.5;

rf = RectPulse(B1, f_rf, 0, tp);


f_rf = gamma*B0;
tp = 1e-3;
alpha=pi/2;
B1 = alpha/(gamma*tp); 


% Tissue parameters
T1GM = 833e-3; % s
T2GM = 83e-3; % s
rhoGM = 0.86;

T1WM = 500e-3; % s
T2WM = 70e-3; % s
rhoWM = 0.77;

FOV = 8e-3;
pixel_size = 1e-3; 

kmaxx = 1/pixel_size; 
kmaxy = 1/pixel_size; 

dkx = 1/FOV;
dky = 1/FOV;

kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid

tau = 4e-3;

N = length(kya)
% Calculate the gradients as the delta_k/(gamma*tau)
Gpexs = kxa/(gamma*tau); 

gx = Gradient([tp tp+tau], {Gpexs 0});


% Frequency encoding gradients
Gfey1 = -2*(kmaxy/(gamma*tau));
Gfey2 = kmaxy/(gamma*tau);

% Need to sample twice times the frequency(Nyqvist)
max_freq = 2/pixel_size

dt = 1/max_freq


%% Test 1 - T1

TE_T1 = -min(T2WM, T2GM) * log(0.93) % 7% signal dephasing
c_w = rhoWM * exp(-TE_T1/T2WM)
c_g = rhoGM * exp(-TE_T1/T2GM)
TR_T1 = abs(log(T1WM*c_g/(T1GM*c_w))*T1GM*T1WM/(T1WM-T1GM))



TE = TE_T1
TR = TR_T1


% same area of both freq encoding pulses
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

sampler = ADC(TE-tau, TE+tau, dt);

[S, ts] = brain_halfmm_pixel(B0, rf, gx, gy, sampler, TR, length(Gpexs));
figure()

[im, x, y] = mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true);


%% Test - T2

TR_T2 = -max(T1WM, T1GM) * log(0.35)

c_w = rhoWM * (1 - exp(-TR_T2/T1WM))
c_g = rhoGM * (1 - exp(-TR_T2/T1GM))
TE_T2 = abs(log(T2WM*c_g/(T2GM*c_w))*T2GM*T2WM/(T1WM-T1GM))



TE = TE_T2
TR = TR_T2


gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);


sampler = ADC(TE-tau, TE+tau, dt);

[S, ts] = brain_halfmm_pixel(B0, rf, gx, gy, sampler, TR, length(Gpexs));
figure()

[im, x, y] = mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true);



%% Test proton weighted

FOV = 8e-3;
pixel_size = 1e-3; 

kmaxx = 1/pixel_size; 
kmaxy = 1/pixel_size; 

dkx = 1/FOV;
dky = 1/FOV;


kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid


TR_P = -max(T1WM, T1GM) * log(0.05)
TE_P = -min(T2WM, T2GM) * log(0.93)


TE = TE_P
TR = TR_P


tau = 4e-3
N = length(kya)
% Calculate the gradients as the delta_k/(gamma*tau)
Gpexs = kxa/(gamma*tau); 

gx = Gradient([tp tp+tau], {Gpexs 0});


% Frequency encoding gradients
Gfey1 = -2*(kmaxy/(gamma*tau));
Gfey2 = kmaxy/(gamma*tau);
% same area of both freq encoding pulses
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

% Need to sample twice times the frequency(Nyqvist)
max_freq = 2/pixel_size

dt = 1/max_freq

sampler = ADC(TE-tau, TE+tau, dt);

[S, ts] = brain_halfmm_pixel(B0, rf, gx, gy, sampler, TR, length(Gpexs));
figure()

[im, x, y] = mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true);

%% Chemical shift


delta = 3.5e-6; % Chemical shift for fat in ppm

kmax = 1 / pixel_size; % Maximum spatial frequency in m^-1


G_fe = kmax / (gamma * tau); 

% Calculate tau for a chemical shift of one pixel
tau = delta * B0 / (G_fe * 1)

%% Test simulation

gammabar = 42.58e6;
gamma = 2*pi*gammabar;
B0 = 1.5;



f_rf = gamma*B0;
tp = 1e-3;
alpha=pi/2;
B1 = alpha/(gamma*tp); 


rf = RectPulse(B1, f_rf, 0, tp);

FOV = 8e-3;
pixel_size = 1e-3; 

iv = waterandfat(3, 1);

tau = 8e-3


TE = 2*tau+tp+0.5*tp; 
TR = 2; 

kmaxx = 1/pixel_size; 
kmaxy = 1/pixel_size; 

dkx = 1/FOV;
dky = 1/FOV;

kxa = -kmaxx:dkx:kmaxx-dkx; % Array
kya = -kmaxy:dky:kmaxy-dky; % Array
[kxg, kyg] = meshgrid(kxa,kya); % Grid

Gpexs = kxa/(gamma*tau); 

gx = Gradient([tp tp+tau], {Gpexs 0});


% Frequency encoding gradients
Gfey1 = -2*(kmaxy/(gamma*tau));
Gfey2 = kmaxy/(gamma*tau);
% same area of both freq encoding pulses
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gfey1 0 Gfey2 0]);

% Need to sample twice times the frequency(Nyqvist)
max_freq = 2/pixel_size;

dt = 1/max_freq;

sampler = ADC(TE-tau, TE+tau, dt);

[S, ts] = seemri(iv, B0, rf, gx, gy, sampler, TR, length(Gpexs), 'PlotKSpace', true);

figure()


[imout, x, y] = mrireconstruct(S, max(kmaxx,kmaxy), 'Plot', true);


