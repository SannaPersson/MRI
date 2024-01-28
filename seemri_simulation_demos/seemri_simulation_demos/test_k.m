%%
N = 8;

x = -N/2:N/2-1;
kx = [-N/2:N/2-1]/N;

%%
k = [0 0 0 0 0 -1i 0 0];
s = fft(ifftshift(k))
plot([real(s); imag(s)]')

