gammabar = 42.58e6;
gamma = 2*pi*gammabar;

B0 = 1;

u = @(x,y,R) double(sqrt(x.^2+y.^2)<=R);

dx = 1;
dy = 1;
x = -3:dx:3;
y = -3:dy:3;
[xs, ys] = meshgrid(x, y);
Mz0 = u(xs, ys, 3);
Mz0(ys>0) = Mz0(ys>0)*0.5;
iv = ImagingVolume(x, y, 0.8, 0.07, Mz0, 'dB0Gamma', 0e-6);


% Sampling
W = 10;
dk_max = 1/W;
kmax = 1/dx/2;
dk = kmax/ceil(kmax/dk_max)

TE = 10e-3;
TR = 5;

% RF pulse
tp = 1e-3;
B1 = pi/2/gamma/tp;
rf = SincPulse([B1 2*B1], B0*gammabar, [0 pi/2], tp, [0 TE/2]);

% Gradients
tau = (TE+tp/2-tp)/2-tp;
ks = -kmax:dk:kmax-dk;
Gpexs = ks/gammabar/tau;
gx = Gradient([tp tp+tau], {Gpexs 0});
ky0 = -kmax;
Gpey = ky0/gammabar/tau
Gfe = kmax/gammabar/tau
gy = Gradient([tp tp+tau TE+tp/2-tau TE+tp/2+tau], [-Gpey 0 Gfe 0]);

% ADC
dt = dk/gammabar/abs(Gfe);
adc = ADC(TE+tp/2-tau, TE+tp/2+tau, dt);

N = length(adc.Times);
actual_kmax = N/2*dk;

dt = dk/gammabar/abs(Gfe);
adc = ADC(TE-tau, TE+tau, dt);

%[S, ts] = seemri(iv, B0, rf, gx, gy, adc, TR, length(Gpexs), ...
%                 'Plot', true, 'TimeStep', dt, ...
%                 'Pause', [0 tp (TE+tp/2)/2 (TE+tp/2)/2+tp]);
[S, ts] = seemri(iv, B0, rf, gx, gy, adc, TR, length(Gpexs), ...
                 'Plot', false);

mrireconstruct(S, kmax, 'Plot', true, 'Filter', 'none');
