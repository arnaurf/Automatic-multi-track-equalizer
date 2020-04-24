clear all; 
%% %%%%%%% INIT VARIABLES
% WINDOW & OVERLAP
winSize = 1024; %window size
overlapR = 0.5; %overlap %
overlapSamples = ceil(overlapR*winSize);
stepL = winSize - overlapSamples;

% LOAD AUDIOS -> loadAudios(directory, winSize, overlapR, t0, tf)
x_original = loadAudios("audios", winSize, overlapR, 60, 65); %t0 and tf optionals

nTracks = size(x_original, 1);
Length = size(x_original.samples{1},2);
fs = x_original.fs{1};

y_filtered = zeros(nTracks, Length); %create empty output

% INIT FILTER
Q = 2;
nF = 10; %number of filters EQ
aF = 10; %number of analisis bands (bins)
fc = 31.25*2.^(0:10-1);
filter = table(zeros(nTracks, nF)+fc, zeros(nTracks, nF), 'VariableNames', ["center", "gain"]);
fc = [1,fc, fs/2-1];

%butterworth fiter at fc(iband)
%b = zeros(length(fc), 5); a = b;
for iband=2:length(fc)-1
    fcs = [fc(iband-1), fc(iband+1)] ./ (fs/2);  
    fcs(2) = min(fc(2),0.9999);
    [b(iband-1,:),a(iband-1,:)] = butter(2,fcs,'bandpass'); 
end

% OTHER VARIABLES
Rt = 3; %max Rank
W = hamming(winSize);
alpha = 0; %alpha value for EMA smoothing

% EMPTY FRAME VARIABLE
frame = table(zeros(nTracks, winSize), zeros(nTracks, nF), zeros(nTracks, nF));
frame.Properties.VariableNames = ["Samples", "MagRes", "Rank"];


%% MAIN LOOP: for each frame of size winSize
nframe = 1;
winPos = 1;
last_t = clock;
tic
while winPos <= Length-winSize
    
    %Read frame for each track
    for iTrack = 1:nTracks
        samples = x_original.samples{iTrack}(winPos:winPos+winSize-1).*W'; %read frame
        
        %%%%%%%%%%%% FEATURES EXTRACTION
        MagRes = getMagRes(samples, b,a);        %get Magnitude of each band
        Rank = getRank(MagRes, aF)';                %get Ranking of most important bands on the frame
        frame(iTrack,:) = table(samples, MagRes, Rank);  %save it on a table
    end
    

%%%%%%%%%%%%%%%%%%%%% MASKING DETECTION
    masking(1).M = zeros(nTracks, nF); %init masking table -> masking(Masker).M(Maskee, band)
    
    for i_masker= 1:nTracks %for each masker
        for i_maskee = 1:nTracks %compare it with each possible maskee
            if i_masker ~= i_maskee %avoid comparing with itself
                for i_band = 1:nF %per each band
                    if (frame.Rank(i_maskee,i_band) <= Rt) && (Rt < frame.Rank(i_masker, i_band)) %Equation (1) on reference 1
                        masking(i_masker).M(i_maskee, i_band) = frame.MagRes(i_masker,i_band) - frame.MagRes(i_maskee,i_band);
                    else
                        masking(i_masker).M(i_maskee, i_band) = 0;
                    end
                end
            else %if i_masker = i_maskee, do:
                masking(i_masker).M(i_maskee, :) = zeros(1, nF);
            end %end if
        end
    end
    
%%%%%%%%%%%%% MASKING SELECTION

    for i_masker = 1:nTracks %for each masker
        
        %Get max masking for the current masker on each band
        Mask = selectMasking(masking(i_masker).M, nF); 
        
        %Smoothing between frames to avoid artefacts
        for i_bin = 1:nF
            if(Mask(i_bin) ~= 0)
                filter.gain(i_masker,i_bin) = EMA(Mask(i_bin), filter.gain(i_masker,i_bin), alpha);
                %For now we have 10 fixed filter centers, no need to smooth
                %filter.center = EMA(center(i_bin), filter.center(i_bin), alpha);
            else
                filter.gain(i_bin) = 0;
                %filter.center(i_bin) = 0;
            end
        end
        
    end
    
%%%%%%%%%%%%%%%% FILTERING
    frame_filtered = zeros(nF, winSize); %init empty filtered frame variable
    for i_track = 1:nTracks
        
        %EQUALIZE   eq_filter(x, fc, Q, winSize, fs)
        frame_filtered = eq_filter(frame.Samples(i_track,:), filter.center(i_track,:), zeros(1, 10)+Q, filter.gain(i_track,:).*(-1), fs);
        
        %Sum the overlapping part of samples with new ones
        y_filtered(i_track,winPos:winPos+winSize-1) = [y_filtered(i_track,winPos:winPos+overlapSamples-1), zeros(1, stepL)] + frame_filtered;
    end

    % First frame alpha is 0.
    if(alpha == 0)
        alpha = exp(-1/(2*fs));
    end
    
    winPos = winPos + stepL;
    nframe =  nframe+1;
    if etime(clock,last_t) >= 1
        fprintf("Frame %i/%i, %.2f%% \n",nframe, ceil(Length/stepL), nframe/ceil(Length/stepL)*100 );
        last_t = clock;
    end
end
t = toc;
fprintf("Elapsed time is %f seconds, performance: %fx\n", t, Length/fs/t);
%% Normalize
% for i = 1:nTracks
%     tracksFilt(i).samples = tracksFilt(i).samples/max(tracksFilt(i).samples)*0.8; %0.8 because want to be sure
% end


%% WRITE AUDIOS
for i = 1:nTracks
    audiowrite("audios_rendered/FILT-"+x_original.fileName{i}, real(y_filtered(i,:)), fs);
end

for i = 1:nTracks %Rewrite original audios bc it could have been cut
    audiowrite("audios_rendered/ORIG-"+x_original.fileName{i}, real(x_original.samples{i,:}), fs);
end

%%
sound(real(y_filtered(1,:)), fs);

%%
clear sound

%% %%%%%%%%%%%%  REFERENCES  %%%%%%%%%%%%%
% [1] Hafezi, S., & Reiss, J. (2015). Autonomous Multitrack Equalization Based on Masking Reduction. Journal of the Audio Engineering Society, 63(5), 312–323. doi: 10.17743/jaes.2015.0021
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%