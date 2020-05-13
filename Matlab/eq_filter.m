function y = eq_filter(x, fc, Q, M, fs)

nozero = find(M);
centers = fc(nozero); %delete 0dB centers
if isempty(nozero)
    y = x;
    return;
end
magnitude = M(nozero);
Qfact = Q(nozero);
N = ones(size(nozero,1),1)'+1;

w0 = centers/(fs/2);
BW = w0/Qfact;
[B,A] = designParamEQ(N,magnitude,w0,BW);
myFilter = dsp.BiquadFilter("SOSMatrixSource","Input port","ScaleValuesInputPort",false);
y = myFilter(x,B,A);
%SOS = [B',[ones(sum(N)/2,1),A']];
%[h, t] = impz(SOS, 48000);
%H = fft(h, size(x,1));
%X = fft(x, size(x,1));
%Y = X.*H;
%y = real(ifft(Y));
%y = filter(B,A,x);
end
