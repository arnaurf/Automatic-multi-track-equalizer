function [mask, mask_curve]  = maskAmount(x1, x2, Fs, plot_on)

    %Parameters
    N = 4096;

    % Other needed variables
    f_kHz = (1:Fs/N:Fs/2)/1000;
    A = 3.64*(f_kHz).^(-0.8) - 6.5*exp(-0.6*(f_kHz - 3.3).^2) + (10^(-3))*(f_kHz).^4;
    %semilogx(0:(Fs/2)/(N/2):Fs/2-1,A,'b'); % Plot in log scale

    % 1kHz signal for SPL normalisation
    Ns = 48000/1000;
    s = sin(2*pi*1000*(0:1/48000:1));
    S = fft(s,Ns);
    S1 = abs(S(1:end/2+1))*2/Ns; %Half of the spectrum (x2 amplitude) and normalize fft (/N)
    fft_max = max(abs(S));

    
    x1w = x1(abs(x1)>0.005);
    x2w = x2;
    
    % Normalize loudness to -10dB (respect the loudest signal)
    x2_loud= integratedLoudness(x2w',Fs);
    dBgain = -10 - x2_loud;
    x1w = x1w*2^(dBgain/6);
    x2w = x2w*2^(dBgain/6);
    %% 
    X1 = fft(x1w(x1w~=0), N);
    X2 = fft(x2w(x2w~=0), N);
    
    %Normalize signal to dB SPL
    fft_spl1 = 96 + 20*log10( abs(X1(1:end/2))*2/N  );
    fft_spl2 = 96 + 20*log10( abs(X2(1:end/2))*2/N  );

    %% Find peaks
    [peak_value, peakPos] = findpeaks(abs(fft_spl2));

    %% 
    big_mask = Schroeder(Fs,N,peakPos'*Fs/N,fft_spl2(peakPos)',14);
    curve = 0;
    for p=1:size(big_mask, 1)
        curve = max(curve,max( A, big_mask(p, :) ) );
    end
    
    curve2 = 0;
    for p=1:size(big_mask, 1)
        curve2 = max(curve2, max( A, big_mask(p, :) ) );
    end

    %%
    f=(0:N/2-1)*Fs/N;
    signal_fRange = fft_spl1>A; %%Only compute masking when signal is audible
    if any(signal_fRange)
        mask_curve = curve2(signal_fRange) - (fft_spl1(signal_fRange));
    else
        mask_curve = 0;
    end
    mask_curve(mask_curve<0) = 0;
    %signal_smr(signal_smr<0) = 0;
    mask = mean(mask_curve);
    mask_curve = {[mask_curve; f(signal_fRange)]};
    %=========PRINT RESULTS=========%
    if plot_on
        figure
        semilogx(f,  fft_spl1, "DisplayName", "Maskee");
        legend('on')
        hold on; grid on;
        semilogx(f,  fft_spl2, "DisplayName", "Masker");
        semilogx(f,  A, "DisplayName", "Threshold in quiet");
        semilogx(f,  curve2, "DisplayName", "Masking Threshold");
        semilogx(f(signal_fRange),  mask_curve, "DisplayName", "Masking amount");
        hold off;
    end
    %===============================%
end