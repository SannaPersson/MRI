%%
tp = 1e-6;

gamma = 2*pi*42.58e6;

B0 = 1;

tp = 1e-6;
B1 = pi/2/gamma/tp;

rf = SincPulse(1*B1, 42.58e6, 0, tp, 0);

iv = ImagingVolume(0, 0, 0.5, 0.2);

[S, ts] = seemri(iv, B0, rf, [], [], ADC(2*tp, tp/100));
