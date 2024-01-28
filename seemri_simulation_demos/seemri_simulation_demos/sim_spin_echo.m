gamma = 2*pi*42.58e6;

B0 = 1;

tp = 1e-4;
B1 = pi/2/gamma/tp;

TE = 0.5e-2;
tau = TE/4;

rf = SincPulse([1 2]*B1, 42.58e6, [0 0], tp, [0 TE/2+tp/4]);

Gx = 0.0004;
gx = Gradient([0 tp 1.5*tp TE/2+tp/4 TE/2+5*tp/4], ...
              [Gx -Gx 0 Gx 0]);
Gy = 0.00001;
gy = Gradient([0 1.6*tp 1.6*tp+tau TE-tau TE+tau], [0 Gy 0 Gy 0]);

iv = ImagingVolume(-5:0.5:5, -5:0.5:5, 1, 0.8, 1, 'dB0Gamma', 1e-6);
%iv = ImagingVolume(-5:1:5, 0, 1, 0.8, 0.1e-6);

[S, ts] = seemri(iv, B0, rf, gx, gy, ADC(TE+tau, 5e-5), 'Layout', 2);

