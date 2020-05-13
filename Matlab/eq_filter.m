function y = eq_filter(x, fc, Q, M, fs)

nozero = find(abs(M)>0.001); %discard small and null filters 
centers = fc(nozero); %delete 0dB centers
if isempty(nozero)
    y = x;
    return;
end
magnitude = M(nozero);
Qfact = Q(nozero);
N = ones(size(nozero,1),2)+1;

%w0 = centers/(fs/2);
%BW = w0./Qfact;
%[B,A] = designParamEQ(N,magnitude,w0,BW);
%myFilter = dsp.BiquadFilter("SOSMatrixSource","Input port","ScaleValuesInputPort",false);
[B, A] = peakFilter(centers(1), Qfact(1), fs, magnitude(1));
for i=2:length(nozero)
    [B2, A2] = peakFilter(centers(i), Qfact(i), fs, magnitude(i));
    B = conv(B,B2); %Cascade filters
    A = conv(A,A2);     
end
y = filter(B,A,x);%myFilter(x,B,A);
%SOS = [B',[ones(sum(N)/2,1),A']];
%[h, t] = impz(SOS, 48000);
%H = fft(h, size(x,1));
%X = fft(x, size(x,1));
%Y = X.*H;
%y = real(ifft(Y));
%y = filter(B,A,x);
end
