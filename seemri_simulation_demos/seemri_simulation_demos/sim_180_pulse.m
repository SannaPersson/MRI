tp = 1e-6;

gamma = 2*pi*42.58e6;

B0 = 1;

tp = 1e-7;
B1 = pi/2/gamma/tp;

TE = 1e-4;
tau = TE/3;

rf = SincPulse(B1, 42.58e6, 0, tp, 0);

Gx = 0.02;
gx = Gradient([0 tp 1.57*tp], [Gx -Gx 0]);

iv = ImagingVolume(-5:1:5, -5:1:5, 0.5, 0.2, 0e-6);
%iv = ImagingVolume(-5:1:5, 0, 0.5, 0.2, 0e-6);

[S, ts] = mri(iv, B0, rf, gx, [], ADC(1.65*tp, 1e-9));
