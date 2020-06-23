function showSpectrum(x, fs, fc, name)

if(size(x, 1)>1)
    x = (x(1,:)+x(2,:))./2;
end

fc2 = fc(2:end-1) + fc(2:end-1)./2;
fc2 = [fc2(1)/2, fc2];

N = fs/fc2(1);
X = abs(fft(x, N));
Xmag = 20*log10(X(1:end/2)*2/N);
plot(linspace(0,fs/2, N/2), movmean( Xmag , 8));

hold on;
stem(fc2, ones(length(fc2),1)*(min(Xmag)), 'b');

set (gca , 'XScale' ,'log' ); grid on;
title(name, 'Interpreter', 'none');
ylim([min(Xmag), 0]);
xticks([0.010,0.020,0.050,0.100,0.200,0.500,1,2,5,10,20].*1000);
xticklabels({'10','20','50','100','200','500','1k','2k','5k','10k','20k'})
xlabel("Frequency (Hz)");
ylabel("Power (dB)");
text(fc(2:end-1), ones(length(fc)-2, 1).*-5, string(0:9))
hold off;
end