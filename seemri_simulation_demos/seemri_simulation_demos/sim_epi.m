gammabar = 42.58e6;
gamma = 2*pi*gammabar;

B0 = 1;

u = @(x,y,R) double(sqrt(x.^2+y.^2)<=R);

dx = 1;
dy = 1;
x = -3:dx:3;
y = -3:dy:3;
[xs, ys] = meshgrid(x, y);
Mz0 = u(xs, ys, 3);
Mz0(ys>0) = Mz0(ys>0)*0.5;
iv = ImagingVolume(x, y, 0.8, 0.05, Mz0);

% RF pulse
tp = 1e-5;
B1 = pi/2/gamma/tp;
rf = SincPulse(B1, B0*gammabar, 0, tp);

% Sampling
W = 8;
dk_max = 1/W;
kmax = 1/dx/2;
dk = kmax/ceil(kmax/dk_max)
kmax = 1/dx/2;
ks = -kmax:dk:kmax-dk;

Ttot = 50e-3;
TE = 30e-3;

Ntot = length(ks)^2;
dt = Ttot/Ntot;
dG = dk/gammabar/dt;


Tcol = (length(ks)-1)*dt
Tstart = TE-dt*(length(ks)^2/2 + length(ks)/2);
Tpe = min(Tstart-tp, 5e-3);

ts = zeros(1, 2*length(ks));
gys = zeros(1, 2*length(ks));
gxs = zeros(1, 2*length(ks));
for ind = 1:length(ks)    
    ts(2*ind-1) = (Tcol+dt)*(ind-1);
    ts(2*ind) = (Tcol+dt)*(ind-1) + Tcol;
    gys(2*ind-1) = (-1)^(ind-1)*dG;
    gys(2*ind) = 0;
    gxs(2*ind-1) = 0;
    if ind < length(ks)
        gxs(2*ind) = dG;
    end
end

gy = Gradient([tp tp+Tpe Tstart+ts], [-kmax/gammabar/Tpe 0 gys]); 
gx = Gradient([tp tp+Tpe Tstart+ts Tstart+ts(end)+dt], ...
              [-kmax/gammabar/Tpe 0 gxs 0]);

adc = ADC(Tstart, Tstart + Ntot*dt, dt);

[S, ts] = seemri(iv, B0, rf, gx, gy, adc, ...
                 'Plot', true, 'TimeStep', dt/4, ...
                 'Pause', [0], ...
                 'PlotKSpace', true);

S2 =reshape(S, 8, 8);    
S2(:,2:2:end) = S2(end:-1:1, 2:2:end);
mrireconstruct(S2, kmax, 'Plot', true)
title('Reconstructed image')
