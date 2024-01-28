%tp = 1e-6;

%gamma = 2*pi*42.58e6;

%B0 = 0.2;

%tp = 5e-7;
%B1 = pi/2/gamma/tp;

iv = ImagingVolume(0, 0, 0.5, 0.2);
rf = RectPulse(0.012, 85.16e5, 0, 5e-7, 0);
[S, ts] = seemri(iv, 0.2, rf, [], [], ADC(7.5e-7, 5e-9), ...
                 'Frame', 'Stationary', 'Pause', 0);
