%gamma = 2*pi*42.58e6;

%tp = 1e-5;
%B1 = pi/2/gamma/tp;

%TE = 0.01;

rf = RectPulse([1 2]*0.00059, 42.58e6, [0 pi/2], 1e-5, [0 0.005]);
iv = ImagingVolume(-4:1:4, -4:1:4, 0.5, 0.02, 1, 'dB0Gamma', 1e-6, ...
                   'PlotScale', 1.6);
[S, ts] = seemri(iv, 1, rf, [], [], ADC(0.02, 0.0001), ...
                 'Pause', [0 1e-5 0.005 0.005+1e-5], 'Layout', 2);

