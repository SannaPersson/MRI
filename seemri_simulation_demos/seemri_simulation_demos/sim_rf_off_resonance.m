%tp = 1e-6;

%gamma = 2*pi*42.58e6;

%B0 = 1;

%tp = 1e-6;
%B1 = pi/2/gamma/tp;

%rf = RectPulse(1*B1, 42.58e6-0.5e6, 0, tp, 0);
rf = RectPulse(0.006, 42.58e6-0.5e6, 0, 1e-6, 0);
iv = ImagingVolume(0, 0, 0.5, 0.2, 1);
[S, ts] = seemri(iv, 1, rf, [], [], ADC(2e-6, 1e-8), 'Pause', 0);
