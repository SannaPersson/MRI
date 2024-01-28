gammabar = 42.58e6;
gamma = 2*pi*gammabar;

iv = ImagingVolume(-2:1, 0, 1e9, 1e9, [2 3 0 1 ], ...
                   'PlotScale', 2);


B0 = 1;
tp = 1e-4;
B1 = pi/2/gamma/tp;
rf = SincPulse(B1, B0*gammabar, 0, tp);

W = 6;
dk_max = 1/W;
dk = dk_max;
kmax = 2;

tau = 1e-3;
Gfe = kmax/gammabar/tau
gx = Gradient([0 tp tp+tau], [0 Gfe 0]);



dt = dk/gammabar/abs(Gfe);
adc = ADC(tp, tp+tau, dt);

[S, ts] = seemri(iv, B0, rf, gx, [], adc, 'TimeStep', dt/10, ...
                 'Pause', [0 tp tp+tau-1.01*dt], 'Layout', 1);

N = length(S);
actual_kmax = N/2*dk;
x2 = (N/2:-1:-N/2+1)/actual_kmax/2;
plot(x2, abs(ifft(fftshift(S').*(-1).^(0:1:N-1))))
title('\fontsize{18}Fourier Transform of Signal')
ylabel('\fontsize{12}|F\{S(t)\}|')

