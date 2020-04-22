function output = getMagRes(x, fcs, fs)

for iband = 2:length(fcs)-1
    %butterworth fiter at fc(iband)
    fc = [fcs(iband-1), fcs(iband+1)] ./ (fs/2);  
    fc(2) = min(fc(2),0.9999);
    [b,a] = butter(2,fc,'bandpass'); 
    %freqz(b,a)
    filtered = filter(b,a,x);

    rms = sqrt(1/size(filtered, 1)*sum(filtered.^2));
    if(rms ~= 0)
        output(iband-1) = 10*log10(rms);
    else
        %output(i) = -100000;
        output(iband-1) = -Inf;
    end
    
end
% 
% function output = getMagRes(x, nAF)
% h = load("butterworth_filter.mat");
% 
% for i = 1:nAF
%     [num, den] = firhalfband(6,  
%     
%     band = conv(x, h.h);
%     rms = sqrt(1/size(band, 1)*sum(band.^2));
%     if(rms ~= 0)
%         output(i) = 10*log10(rms);
%     else
%         %output(i) = -100000;
%         output(i) = -Inf;
% end
%     
% end