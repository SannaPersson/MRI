
gamma = 2*pi*42.58e6;

[S, ts] = seemri(ImagingVolume(-5:1:5, -5:1:5, 0.1, 0.02), 1, ...
              SincPulse(pi/2/gamma/1e-6, 42.58e6-100, 0, 1e-6), ...
              Gradient([0 1e-6 1.51e-6], [0.05 0 0]), [], ...
              ADC(1.51e-6, 1e-8));
