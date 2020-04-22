function output = getMagRes(x, nAF)
h = load("PB_filters.mat");

for i = 1:nAF
    band = conv(x, h.imp{i});
    rms = sqrt(1/size(band, 1)*sum(band.^2));
    if(rms ~= 0)
        output(i) = 10*log10(rms);
    else
        %output(i) = -100000;
        output(i) = -Inf;
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