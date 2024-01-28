gamma = 2*pi*42.58e6;

B0 = 1;

tp = 1e-6;
B1 = pi/2/gamma/tp;

rf = SincPulse(B1, B0*42.58e6, 0, tp, 0);

tissues = [2 3 4 5];
tissuenames = {'Bg', 'CSF', 'Grey matter', 'White matter', 'Fat'};

tissue2T1 = [0 2569 833 500 350 0 569 0 0 500 2569 500]*1e-3;
tissue2T2 = [0 329 83 70 70 47 329 0 0 70 329 70]*1e-3;
tissue2PD = [0 1 0.86 0.77 1 1 1 0 0 0.77 1 0.77];

iv = ImagingVolume(1:length(tissues), 0, ...
                   tissue2T1(tissues), tissue2T2(tissues), ...
                   tissue2PD(tissues), ...
                   'Annotation', ...
                   @(n,T1,T2,Mz0) sprintf(...
                       '%s\nT1 = %g s\n T2 = %g s\nPD = %g', ...
                       tissuenames{tissues(n)}, T1, T2,Mz0));

[S, ts] = seemri(iv, B0, rf, [], [], ADC(0.5, 0.002), ...
                 'Pause', 0);
