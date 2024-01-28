%gamma = 2*pi*42.58e6;

%tp = 1e-5;
%B1 = pi/2/gamma/tp;

%TE = 0.01;

rf = RectPulse(0.00059, 42.58e6, 0, 1e-5, 0);
iv = ImagingVolume(-4:1:4, -4:1:4, 0.5, 0.02, 1, 'dB0Gamma', 1e-6, ...
                   'PlotScale', 1.75);
[S, ts] = seemri(iv, 1, rf, [], [], ADC(0.02, 0.00015), 'Pause', 0, ...
                 'Layout', 2);

