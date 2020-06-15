%clear all;
function signal_smr = smr(x, Fs)
    % Audio reading
    %[x, Fs] = audioread("drum_stem.wav");
    if size(x, 2)==2
        x = sum(x, 2);
    end

    %Parameters
    N = 1024;

    % Other needed variables
    f_kHz = (1:Fs/N:Fs/2)/1000;
    A = 3.64*(f_kHz).^(-0.8) - 6.5*exp(-0.6*(f_kHz - 3.3).^2) + (10^(-3))*(f_kHz).^4;
    %semilogx(0:(Fs/2)/(N/2):Fs/2-1,A,'b'); % Plot in log scale

    % 1kHz signal for SPL normalisation
    s = sin(2*pi*1000*0:10);
    S = fft(s);
    fft_max = max(abs(S));

    %% 
    X = fft(x, N);

    %Normalize signal to dB SPL
    fft_spl = 96 + 20*log10( abs(X(1:end/2))./ fft_max );

    %% Find peaks
    [peak_value, peakPos] = findpeaks(abs(fft_spl));

    %% 
    big_mask = Schroeder(Fs,N,peakPos*Fs/N,fft_spl(peakPos),14);
    curve = 0;
    for p=1:size(big_mask, 1)
        curve = max(curve,max( A, big_mask(p, :) ) );
    end

    %%
    signal_smr = fft_spl' - curve;
    signal_smr(signal_smr<0) = 0;
end