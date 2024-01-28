gamma = 2*pi*42.58e6;

B0 = 1;

tp = 1e-6;
B1 = pi/2/gamma/tp;

rf = SincPulse(B1, B0*42.58e6, 0, tp, 0);

iv = ImagingVolume(-1:1:1, 0, 0.5, [0.05 0.1 0.5], 1, ...
                   'Annotation', @(n,T1,T2,PD) sprintf('T2 = %g s', T2));

[S, ts] = seemri(iv, B0, rf, [], [], ADC(0.2, 0.002), ...
                 'Pause', 0, 'Layout', 2);
