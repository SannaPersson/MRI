gammabar = 42.58e6;
gamma = 2*pi*gammabar;

W = 10;
dk_max = 1/W;
dk = dk_max;

kmax = 1;

TE = 1e-3;
tau = TE/4;

kx0 = 0;
ky0 = -kmax;

Gpex = kx0/gammabar/tau
Gpey = ky0/gammabar/tau
Gfe = kmax/gammabar/tau

dt = dk/gammabar/abs(Gfe);

B0 = 1;

tp = 1e-5;
B1 = pi/2/gamma/tp;

rf = SincPulse(B1, B0*gammabar, -pi/2, tp);
gx = Gradient([0 tp tp+tau], [0 Gpex 0]);
gy = Gradient([0 tp tp+tau TE-tau TE+tau], [0 Gpey 0 Gfe 0]);

x = -4:0.5:4;
y = -4:0.5:4;
[xs, ys] = meshgrid(x, y);
rs = sqrt(xs.^2+ys.^2);
Mz0 = double(rs <= 4)*(x(2)-x(1))*(y(2)-y(1));
Mz0(ys>0) = Mz0(ys>0)*0.5;

iv = ImagingVolume(x, y, 1, 0.8, Mz0);

adc = ADC(TE-tau, TE+tau, dt);

N = length(adc.Times);
actual_kmax = N/2*dk;

% 2D imaging
TR = 2;
S2 = zeros(length(adc.Times));
iv = ImagingVolume(x, y, 1, 0.8, Mz0);
ind = 1;
for kx0 = -kmax:dk:kmax-dk
    fprintf(num2str(ind, '%d,'));
    ky0 = -kmax;
    
    Gpex = kx0/gammabar/tau;
    Gpey = ky0/gammabar/tau;

    Gpexs(ind) = Gpex;
    
    gx = Gradient2([0 tp tp+tau], [0 Gpex 0]);
    gy = Gradient2([0 tp tp+tau TE-tau TE+tau], [0 Gpey 0 Gfe 0]);
    
    [S, ts] = seemri(iv, B0, rf, gx, gy, adc, TR, 'Plot', true);
    %seemri(iv, B0, [], [], [], ADC(TR-TE+tau,TR-TE+tau), 'Plot', false);
    
    S2(:,ind) = S;
    ind = ind + 1;
end
fprintf('\n');

figure(3)
[M,N] = size(S2);
x2 = (-N/2:1:N/2-1)/actual_kmax/2;
y2 = (-M/2:1:M/2-1)/actual_kmax/2;
%gk = gausskern2(3.17);
%S3 = gk(1:end-1,1:end-1).*S2;
imagesc(x2, y2, abs(ifft2(fftshift(S2...
                                   .*((-1).^((1:M)'*ones(1,N) ...
                                             +ones(M,1)*(1:N)))))))




W = 10;
dk_max = 1/W;
dk = dk_max;

Gpexs = (-kmax:dk:kmax-dk)/gammabar/tau;
gx = Gradient2([tp tp+tau], {Gpexs 0});
ky0 = -kmax;
Gpey = ky0/gammabar/tau
Gfe = kmax/gammabar/tau
gy = Gradient2([tp tp+tau TE-tau TE+tau], [Gpey 0 Gfe 0]);

dt = dk/gammabar/abs(Gfe);
adc = ADC(TE-tau, TE+tau, dt);

iv.toEquilibrium();
[S, ts] = seemri(iv, B0, rf, gx, gy, adc, TR, length(Gpexs), 'Plot', true);

figure
[M,N] = size(S);
x2 = (-N/2:1:N/2-1)/actual_kmax/2;
y2 = (-M/2:1:M/2-1)/actual_kmax/2;
%gk = gausskern2(3.17);
%S3 = gk(1:end-1,1:end-1).*S2;
imagesc(x2, y2, abs(ifft2(fftshift(S...
                                   .*((-1).^((1:M)'*ones(1,N) ...
                                             +ones(M,1)*(1:N)))))))
