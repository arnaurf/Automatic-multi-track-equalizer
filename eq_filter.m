function y = eq_filter(x, fc, Q, M, fs)

N = [6,6,6,6,6,6,6,6,6,6];
[B,A] = designParamEQ(N,M,fc/(fs/2),fc/(fs/2)./Q, 'sos');
SOS = [B',[ones(sum(N)/2,1),A']];

[h, t] = impz(SOS, 48000);
H = fft(h, size(x,1));
X = fft(x, size(x,1));
Y = X.*H;
y = real(ifft(Y));

end