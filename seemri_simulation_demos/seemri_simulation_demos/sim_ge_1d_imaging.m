gammabar = 42.58e6;
gamma = 2*pi*gammabar;

W = 10;
dk_max = 1/W;
dk = dk_max;

kmax = 2;

TE = 1e-3;
tau = TE/4;

kx0 = 0;
ky0 = -kmax;

Gpex = kx0/gammabar/tau/2
Gpey = ky0/gammabar/tau/2
Gfe = kmax/gammabar/tau/2

dt = dk/gammabar/abs(Gfe);

B0 = 1;

tp = 1e-5;
B1 = pi/2/gamma/tp;

rf = SincPulse(B1, B0*gammabar, -pi/2, tp);
gx = Gradient([0 tp tp+tau], [0 Gpex 0]);
gy = Gradient([0 tp tp+tau TE-tau TE+tau], [0 Gpey 0 Gfe 0]);

x = -4:0.5:4;
y = -4:0.5:4;
[xs, ys] = meshgrid(x, y);
rs = sqrt(xs.^2+ys.^2);
Mz0 = double(rs <= 4)*(x(2)-x(1))*(y(2)-y(1));
Mz0(ys>0) = Mz0(ys>0)*0.5;

iv = ImagingVolume(x, y, 1, 0.8, Mz0);

adc = ADC(TE-tau, TE+tau, dt);

figure(1)
[S, ts] = seemri(iv, B0, rf, gx, gy, adc, 'TimeStep', dt);
S = S';

figure(2)
N = length(S);
actual_kmax = N/2*dk;
x2 = (-N/2:1:N/2-1)/actual_kmax/2;
plot(x2, abs(ifft(fftshift(S).*(-1).^(0:1:N-1))), ...
     x, sum(reshape(iv.Mz0', length(x), length(y)), 2))
legend('Reconstructed', 'Real')

