gamma = 2*pi*42.58e6;

B0 = 1;

tp = 10e-6;
B1 = pi/2/gamma/tp

B1env = @(t) B1.*(t >= 0).*(t <= tp);

rf = RFPulse(B1env, B0*42.58e6-1e2, 0, tp);

[S, ts] = run_bloch(B0, rf, 0.1, 0.0002);

