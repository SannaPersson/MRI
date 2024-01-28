gammabar = 42.58e6;
gamma = 2*pi*gammabar;

iv = ImagingVolume(-2:1, 0, 0.1, 0.05, [2 3 0 1 ], ...
                   'PlotScale', 2);


B0 = 1;
tp = 1e-4;
B1 = pi/2/gamma/tp;
rf = SincPulse(B1, B0*gammabar, 0, tp);

W = 5;
dk_max = 1/W;
dk = dk_max;
kmax = 1;

tau = 1e-3;
Gpes = (0:dk:kmax-dk)/gammabar/tau
gx = Gradient([0 tp tp+tau], {0 Gpes 0});

dt = dk/gammabar/abs(Gfe);
adc = ADC(tp+tau, tp+2*tau, dt);

figure(1)
[S, ts] = seemri(iv, B0, rf, gx, [], adc, 1, length(Gpes), ...
                 'TimeStep', dt, ...
                 'Pause', [0 tp tp+tau tp+2*tau tp+2*tau-1.01*dt], ...
                 'Layout', 1);
%, ...
%                 'Plot', false);

N = length(S);
actual_kmax = kmax;

figure(2)
subplot(1,2,1)
plot(ts, abs(S))
xlabel('\fontsize{18}Time (s)');
ylabel('\fontsize{18}|S(t)|');
title('\fontsize{18}Phase-Encoded Signals')

subplot(1,2,2)


s = S(1,:)*exp(i*pi/2);
%plot(abs(ifft(fftshift(S(1,:).*(-1).^(0:size(S,2)-1)))))             
x2 = (-N/2:1:N/2-1)/actual_kmax/2;
plot(x2, abs(ifft(fftshift([conj(s(1)) conj(s(end:-1:2)) s(1:end)].*(-1).^(0:2*size(S,2)-1)))))
title('\fontsize{18}Fourier Transform of Combined Signals')


% X2 = (-N/2:1:N/2-1)/actual_kmax/2;
% plot(x2, abs(ifft(fftshift(S(1,:)).*(-1).^(0:1:N-1))))
% title('\fontsize{18}Fourier Transform of Combined Signals')

% plot(abs(ifft(fftshift(S(1,:).*(-1).^(1:size(S,2))))))
