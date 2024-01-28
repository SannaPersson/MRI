

%[ph, xph, yph, zph] = brainwebphantom(4);
%slice = 170;

[current_slice, xph, yph] = ph4_170();

subsamp = 1;
x = xph(1:subsamp:end);
y = yph(1:subsamp:end);

% For subject >= 4: (0=Background, 1=CSF, 2=Gray Matter, 3=White
% Matter, 4=Fat, 5=Muscle, 6=Muscle/Skin, 7=Skull, 8=vessels, 9=around
% fat, 10 =dura matter, 11=bone marrow)
% step 0.5 mm
tissue2T1 = [0 2569 833 500 350 0 569 0 0 500 2569 500]*1e-3;
tissue2T2 = [0 329 83 70 70 47 329 0 0 70 329 70]*1e-3;
tissue2PD = [0 1 0.86 0.77 1 1 1 0 0 0.77 1 0.77];

T1 = tissue2T1(current_slice(1:subsamp:end,1:subsamp:end)+1);
T2 = tissue2T2(current_slice(1:subsamp:end,1:subsamp:end)+1);
PD = tissue2PD(current_slice(1:subsamp:end,1:subsamp:end)+1);
dB0 = -3.35e-6*(current_slice(1:subsamp:end,1:subsamp:end) == 4);

iv = ImagingVolume(x, y, T1, T2, PD, 'dB0', dB0);


gammabar = 42.58e6;
gamma = 2*pi*gammabar;


B0 = 1;

% RF pulse
tp = 1e-5;
B1 = pi/2/gamma/tp;
rf = RectPulse(B1, B0*gammabar, -pi/2, tp);


% Sampling
W = 220;
dk_max = 1/W;
kmax = 1/(0.5*subsamp)/2
dk = kmax/ceil(kmax/dk_max)

% Gradients
TE = 80e-3;
TR = 5;
tau = (TE-tp)/2/10;
ks = -kmax:dk:kmax-dk;
Gpexs = ks/gammabar/tau;
gx = Gradient([tp tp+tau], {Gpexs 0});
ky0 = -kmax;
Gpey = ky0/gammabar/tau
Gfe = kmax/gammabar/tau
gy = Gradient([tp tp+tau TE-tau TE+tau], [Gpey 0 Gfe 0]);

% ADC
dt = dk/gammabar/abs(Gpey);
adc = ADC(TE-tau, TE+tau, dt);

N = length(adc.Times);
actual_kmax = N/2*dk;

dt = dk/gammabar/abs(Gfe);
adc = ADC(TE-tau, TE+tau, dt);

TR = 10;
[S, ts] = seemri(iv, B0, rf, gx, gy, adc, TR, length(Gpexs), 'Plot', false);

[im, x2, y2] = mrireconstruct(S, actual_kmax, 'Plot', true);
%imagesc(x2, y2, im)
%colormap gray
%axis image
