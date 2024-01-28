tp = 1e-6;

gamma = 2*pi*42.58e6;

B0 = 1;

tp = 1e-6;
B1 = pi/2/gamma/tp;

tau = 5e-6;

rf = SincPulse(1*B1, 42.58e6, 0, tp, 0);


Gx = 0.01;
gx = Gradient([0 tp tau], [0 Gx 0]);

iv = ImagingVolume(-4:1:4, -4:1:4, 0.5, 0.2, 100e-6);

[S, ts] = seemri(iv, B0, rf, gx, [], ADC(tau, tau/200), ...
                 'Pause', [0 tp], 'Layout', 2);
