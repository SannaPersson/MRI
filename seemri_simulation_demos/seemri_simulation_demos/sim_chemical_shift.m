gammabar = 42.58e6;
gamma = 2*pi*gammabar;


B0 = 1;

x = -4:0.25:4;
y = -4:0.25:4;
[xs, ys] = meshgrid(x, y);
rs = sqrt(xs.^2+ys.^2);
Mz0 = double(rs <= 4);
Mz0(ys>0) = Mz0(ys>0)*0.5;
iv = ImagingVolume(x, y, 1, 0.8, Mz0, 'dB0', -(Mz0==1)*3.35e-6);

tp = 1e-5;
B1 = pi/2/gamma/tp;
rf = RectPulse(B1, B0*gammabar, -pi/2, tp);

% Sampling
W = 10;
dk_max = 1/W;
dk = dk_max;
kmax = 2;

TEs = 0.005:0.005:0.015;
for ind = 1:length(TEs)
% Gradients
TE = TEs(ind);
TR = 5;
tau = (TE-tp)/2;
ks = -kmax:dk:kmax-dk*1e-3;
[kxs, kys] = meshgrid(ks, ks);
Gpexs = ks/gammabar/tau;
gx = Gradient([tp tp+tau], {Gpexs 0});
ky0 = -kmax;
Gpey = ky0/gammabar/tau
Gfe = kmax/gammabar/tau
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gpey 0 Gfe 0]);

% ADC
dt = dk/gammabar/abs(Gfe);
adc = ADC(TE-tau, TE+tau, dt);

N = length(adc.Times);
actual_kmax = N/2*dk;

dt = dk/gammabar/abs(Gfe);
adc = ADC(TE-tau, TE+tau, dt);

TR = 5;
[S, ts] = seemri(iv, B0, rf, gy, gx, adc, TR, length(Gpexs), 'Plot', false);

[M,N] = size(S);
x2 = (-N/2:1:N/2-1)/actual_kmax/2;
y2 = (-M/2:1:M/2-1)/actual_kmax/2;
%gk = gausskern2(3.17);
%S2 = gk(1:end-1,1:end-1).*S;

% Hamming window
Hk = 0.54 + 0.46*cos(pi*(ks)/kmax);
subplot(1, length(TEs), ind)
imagesc(x2, y2, abs(ifft2(fftshift((S.*(Hk'*Hk)) ...
                                   .*((-1).^((1:M)'*ones(1,N) ...
                                             +ones(M,1)*(1:N)))))))
axis image
title(sprintf('TE = %g s', TE));
colormap gray
end
