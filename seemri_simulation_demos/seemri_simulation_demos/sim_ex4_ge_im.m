gammabar = 42.58e6;
gamma = 2*pi*gammabar;

B0 = 1.5;

% Imaging volume
iv = disc(3,1);

% Imaging parameters
TE = 10e-3;
TR = 2;
tp = 1e-3;
tau = 4e-3;

% RF pulse
B1 = pi/2/gamma/tp;
rf = RectPulse(B1, B0*gammabar, 0, tp);

% k-space sampling
kmax = 1/2/1;
W = 8;
dk_max = 1/W;
dk = kmax/ceil(kmax/dk_max)
ks = -kmax:dk:kmax-dk;

% Gradients
Gpes = ks/gammabar/tau;
Gfe = kmax/gammabar/tau
gx = Gradient([tp tp+tau], {Gpes 0});
gy = Gradient([tp tp+tau TE-tau TE+tau], [-Gfe 0 Gfe 0]);

% ADC
dt = dk/gammabar/abs(Gfe);
adc = ADC(TE-tau, TE+tau, dt);

[S, ts] = seemri(iv, B0, rf, gx, gy, adc, TR, length(Gpes), ...
                 'TimeStep', dt/4, 'PlotKSpace', true);

mrireconstruct(S, kmax, 'Plot', true, 'Filter', 'none');
