function output = getMagRes(x, b, a)

for iband = 1:size(b, 1)
    filtered = filter(b(iband,:),a(iband,:),x);
    
    rms = sqrt(1/size(filtered, 2)*sum(filtered.^2));
    if(rms ~= 0)
        output(iband) = 20*log10(rms);
    else
        output(iband) = -Inf;
    end
    
end

end
