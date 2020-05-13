clear all

%a = [1,0.2];
%b = [1];
x = 1:10;

%x = sin(0.01*0:399);
%y = filter(b, a, x)
x = (0:399)*0.01;
fs = 48000;
Q = 2;
nF = 10; %number of filters EQ
aF = 10; %number of analisis bands (bins)
fc = 31.25*2.^(0:10-1);
fc = [1,fc, fs/2-1];

%butterworth fiter at fc(iband)
%b = zeros(length(fc), 5); a = b;
for iband=2:length(fc)-1
    fcs = [fc(iband-1), fc(iband+1)] ./ (fs/2);  
    fcs(2) = min(fcs(2),0.9999);
    [b(iband-1,:),a(iband-1,:)] = butter(2,fcs,'bandpass'); 
end

a = a(1,:);
b = b(1,:);
for i = 1:size(b,2)
    b(i) = b(i) / a(1);
    a(i) = a(i) / a(1);
end

Y = zeros(length(x),1);
    
for n = 1:length(x)
    auxX = 0;
    auxY = 0;
    for m = 0:size(b,2)-1
        if(n-m > 0)
            auxX = auxX + x(n - m) * b(m+1);
        end
    end
    for m=1:size(b,2)-1
        if (n - m >= 1)
            auxY = auxY + Y(n - m) * a(m+1);
        end
    end
    Y(n) = (auxX - auxY) / a(1);
 end