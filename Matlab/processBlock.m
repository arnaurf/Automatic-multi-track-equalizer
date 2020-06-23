clear all; 
%% %%%%%%% INIT VARIABLES

S = 4;
Q = 5;
nF = 5; %number of filters EQ
alpha_user = 0.99;


% WINDOW & OVERLAP
winSize = 1024; %window size
overlapR = 0.5; %overlap %
overlapSamples = ceil(overlapR*winSize);
stepL = winSize - overlapSamples;

% LOAD AUDIOS -> loadAudios(directory, winSize, overlapR, t0, tf)
x_original = loadAudios("audios", winSize, overlapR); %t0 and tf optionals

nTracks = size(x_original, 1);
Length = size(x_original.samples{1},2);
fs = x_original.fs{1};


y_filtered = cell(nTracks,1);
for i=1:nTracks
    y_filtered{i} = zeros(size(x_original.samples{i})); %create empty output
end
        
% INIT FILTER

aF = 10; %number of analisis bands (bins)
fc = 31.25*2.^(0:10-1);
filter = table(zeros(nTracks, aF)+fc, zeros(nTracks, aF), 'VariableNames', ["center", "gain"]);
fc = [1,fc, fs/2-1];

%butterworth fiter at fc(iband)
%b = zeros(length(fc), 5); a = b;
for iband=2:length(fc)-1
    fcs = [fc(iband-1), fc(iband+1)] ./ (fs/2);  
    fcs(2) = min(fcs(2),0.9999);
    [b(iband-1,:),a(iband-1,:)] = butter(2,fcs,'bandpass'); 
end

% OTHER VARIABLES
Rt = 3; %max Rank
W = (hann(winSize));
% W = hamming(winSize);
% W(winSize) = 0;


alpha = 0; %alpha value for EMA smoothing


% EMPTY FRAME VARIABLE
frame = table(zeros(nTracks, winSize), zeros(nTracks, aF), zeros(nTracks, aF));
frame.Properties.VariableNames = ["Samples", "MagRes", "Rank"];

frameSamples = cell(nTracks,1);
frameMag = zeros(nTracks, aF);
frameRank = zeros(nTracks, aF);

%% MAIN LOOP: for each frame of size winSize
nframe = 1;
winPos = 1;
last_t = clock;
tic
while winPos <= Length-winSize
    
    %Read frame for each track
    for iTrack = 1:nTracks
        frameSamples{iTrack} = x_original.samples{iTrack}(:,winPos:winPos+winSize-1); %read frame
        samples_mono = sum(frameSamples{iTrack}, 1)/2;
        
        if any(samples_mono) && rms(samples_mono) >0.003 %Noise gate
            %%%%%%%%%%%% FEATURES EXTRACTION
            frameMag(iTrack,:) = getMagRes(samples_mono, b, a);        %get Magnitude of each band (dBSPL)
            frameRank(iTrack,:) = getRank(frameMag(iTrack,:), aF)';                %get Ranking of most important bands on the frame
            %frame(iTrack,:) = cell(samples, MagRes, Rank);  %save it on a table
        else
            frameMag(iTrack,:) = ones(1, aF)*(-Inf);
            frameRank(iTrack,:) = zeros(1,aF); %If all bands are marked as essential no masking will be false detected
        end
    end
    

%%%%%%%%%%%%%%%%%%%%% MASKING DETECTION

    M = cell(nTracks);
    for i_masker= 1:nTracks %for each masker
        for i_maskee = 1:nTracks %compare it with each possible maskee
            if i_masker ~= i_maskee %avoid comparing with itself
                for i_band = 1:aF %per each band
                    if (frameRank(i_maskee,i_band) <= Rt) && (Rt < frameRank(i_masker, i_band)) ...
                            && frameMag(i_masker,i_band) ~= -Inf && frameMag(i_maskee,i_band) ~= -Inf %Equation (1) on reference [1]
                        M{i_masker, i_maskee}(i_band) = frameMag(i_masker,i_band) - frameMag(i_maskee,i_band);
                    else
                        M{i_masker, i_maskee}(i_band) = 0;
                    end
                end
            else %if i_masker = i_maskee, do:
                 M{i_masker, i_maskee} = zeros(1, aF);
            end %end if
        end
    end
    
%%%%%%%%%%%%% MASKING SELECTION

    for i_masker = 1:nTracks %for each masker
        
        %If the track is disabled (noise gate) just don't change the
        %filter so when it comes back it has same filtering
        if any(frameRank(i_masker,:))%Check if the track is enabled (noise gate)
            %Get max masking for the current masker on each band
            Mask = selectMasking(cell2mat(M(i_masker,:)')*(-1), aF); 

            new_mask = Mask;
            [val, ind] = sort(new_mask);
            new_mask(ind(1:end-nF)) = 0;

            %Smoothing between frames to avoid artefacts
            for i_bin = 1:aF
                if(Mask(i_bin) ~= 0)
                    filter.gain(i_masker,i_bin) = EMA(new_mask(i_bin)*(-1), filter.gain(i_masker,i_bin), alpha);
                    %For now we have 10 fixed filter centers, no need to smooth
                    %filter.center = EMA(center(i_bin), filter.center(i_bin), alpha);     
                else 
                    filter.gain(i_masker,i_bin) = EMA(0, filter.gain(i_masker,i_bin), alpha);
                end

            end
        end
    end
    
%%%%%%%%%%%%%%%% FILTERING
    for i_track = 1:nTracks
        gain = filter.gain(i_track,:).*S;
        [val, ind] = sort(gain);
        gain(ind(nF+1:end)) = 0;
        %EQUALIZE   eq_filter(x, fc, Q, winSize, fs)
        frame_filtered = eq_filter(frameSamples{i_track}, filter.center(i_track,:), zeros(1, 10)+Q, gain, fs);
        
        %Sum the overlapping part of samples with new ones
        y_filtered{i_track}(:,winPos:winPos+winSize-1) = [y_filtered{i_track}(:,winPos:winPos+overlapSamples-1), zeros(size(frame_filtered,1),stepL)] + frame_filtered.*W';
    end

    % First frame alpha is 0.
    if(alpha == 0)
        alpha = alpha_user;%exp(-1/(2*fs));
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


%% WRITE AUDIOS
y_mixdown = 0;
for i = 1:nTracks
    audiowrite("audios_rendered/FILT-"+x_original.fileName{i}, real(y_filtered{i})', fs);
    y_mixdown = y_mixdown + y_filtered{i};
end

%Original signal
x_mixdown = 0;
for i = 1:nTracks %Rewrite original audios bc it could have been cut
    audiowrite("audios_rendered/ORIG-"+x_original.fileName{i}, real(x_original.samples{i,:}'), fs);
   x_mixdown = x_mixdown + x_original.samples{i};
end

% Equal the loudness
y_loud= integratedLoudness(sum(y_mixdown,1)',fs);
dBgain = integratedLoudness(sum(x_mixdown,1)',fs) - y_loud;
y_mixdown = y_mixdown*2^(dBgain/6);

%Normalise respect max value of both audios
max_value = max( max(y_mixdown(:)), max(x_mixdown(:)));
y_mixdown = y_mixdown./max_value * 0.95;
x_mixdown = x_mixdown./max_value * 0.95;

%Write mixdown version and normalized
audiowrite("audios_rendered/mixdownORIG-"+x_original.fileName{i}, real(x_mixdown)', fs);
audiowrite("audios_rendered/mixdownFILT-"+x_original.fileName{i}, real(y_mixdown)', fs);


%% 
t = 1:nTracks;
mask_amount = zeros(nTracks, 2); %nTracks raws, and two columns: before and after the algorithm
mask_curves = cell(nTracks, 2);
for iT=1:nTracks
    % Mono Track iT, raw audios
    x1 = sum(x_original.samples{iT}, 1)/size(x_original.samples{iT},1);        

    % To mono: all the other tracks
    x2 = zeros(1,length(x1));
    for iM=t(t~=iT) %mono conversion
        x2 = x2 + sum(x_original.samples{iM}, 1)/size(x_original.samples{iM},1);
    end
    
    % -----------------------------------------------------------
    
    %Compute mask_amount for that track vs the others
    [mask_amount(iT,1), mask_curves(iT,1)] = maskAmount(x1, x2, fs, false);

    % Same code but with the filtered versions
    x1 = sum(y_filtered{iT}, 1)/size(y_filtered{iT},1);        %Track iT
    
    x2 = zeros(1,length(x1));
    for iM=t(t~=iT) %mono conversion
        x2 = x2 + sum(y_filtered{iM}, 1)/size(y_filtered{iM},1);
    end
    
    %Compute mask_amount for that track vs the others
    [mask_amount(iT,2), mask_curves(iT,2)] = maskAmount(x1, x2, fs, false);
   
end

improv_track = mask_amount(:,2)./mask_amount(:,1)
improv_overall = mean(improv_track )

%% %%%%%%%%%%%%  REFERENCES  %%%%%%%%%%%%%
% [1] Hafezi, S., & Reiss, J. (2015). Autonomous Multitrack Equalization Based on Masking Reduction. Journal of the Audio Engineering Society, 63(5), 312–323. doi: 10.17743/jaes.2015.0021
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
