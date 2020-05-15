function y = eq_filter(x, fc, Q, M, fs)

nozero = find(abs(M)>0.001); %discard small and null filters 

if isempty(nozero) %if there are any filter with magnitude>0.001, dont process anything
    y = x;
    return;
end

%Select non-zero centers, magnitudes and Q.
centers = fc(nozero); 
magnitude = M(nozero);
Qfact = Q(nozero);

%Cascade filters
[B, A] = peakFilter(centers(1), Qfact(1), fs, magnitude(1));
y = filter(B,A,x);
for i=2:length(nozero)
    [B2, A2] = peakFilter(centers(i), Qfact(i), fs, magnitude(i));
    y = filter(B2,A2,y);  
end
%%%%%%%% Visualize filter
%visualizeFilter(B,A,fs);
%%%%%%%%%%%%%%%%%%%%%%%%%

end


function visualizeFilter(B,A, fs)
    h = fvtool(B,A);
    h.Fs = fs;
    h.FrequencyRange='Specify freq. vector';
    h.FrequencyVector = linspace(10,fs/2,8192);
    h.FrequencyScale = 'Log';
    xticks([0.010,0.020,0.050,0.100,0.200,0.500,1,2,5,10,20]) %Ticks are actually in kHz!
    xticklabels({'10','20','50','100','200','500','1k','2k','5k','10k','20k'})
    xlabel("Frequency (Hz)")
end