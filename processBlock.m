clear all; 
tic

%test git push
%testing editing from gitHub web interface
%Init variables
nF = 10; %number of filters EQ

for n = 1:nF
    filter_d(n).M = 0; %+0dB iniciales
    filter_d(n).center = 31.25*2.^(n-1); %31.25Hz, 62.5Hz, 125Hz etc...
end

Q = 2;
alpha = 0;
aF = 10; %number of analisis bands (bins)
M = 1024; %window size
Rt = 3; %max Rank
overlapR = 0.5; %overlap %

%%%%%%%%%%%% FEATURE EXTRACTION
%Load audios

files = dir("audios");
j = 1;
for i = 1:size(files, 1) %reading all files on "audios" folder
    if files(i).isdir == 0
        [tracks(j).samples, tracks(j).fs] = audioread("audios/"+string(files(i).name));
        tracks(j).samples = tracks(j).samples(9500000:10000000, :); %read only a section of the audio
        tracks(j).samples(end:ceil(size(tracks(j).samples,1)/M)*M+M*overlapR, :) = 0; %Fill input with 0 to fit all windows
        
        tracks(j).name = files(i).name;
        
        tracksFilt(j).samples = zeros(size(tracks(j).samples,1)+M*overlapR,1); %create empty output
        
        if(size(tracks(j).samples, 2) > 1) %Stereo to mono
            tracks(j).samples = sum(tracks(j).samples, 2) / size(tracks(j).samples, 2);
        end
        j = j + 1;
    end    
end

nTracks = size(tracks, 2);
fs = tracks(1).fs;

%for each frame of size M
W = hamming(M);
frame = 1;
w_pos = 1;

while w_pos <= size(tracks(1).samples, 1)-M*overlapR
    
    %Read frame for each track
    for i = 1:nTracks
        tracksFrame(i).samples = tracks(i).samples(w_pos:w_pos+M-1).*W; 
        
        %%%%%%%%%%%% FEATURES EXTRACTION
        tracksFrame(i).MagRes = getMagRes(tracksFrame(i).samples, aF);
        tracksFrame(i).Rank = getRank(tracksFrame(i).MagRes, aF);
    end
    
%%

%%%%%%%%%%%%%%%%%%%%% MASKING DETECTION
    masking = zeros(nTracks, nTracks, nF);
    for masker = 1:nTracks
        for maskee = 1:nTracks
            if masker ~= maskee
                
                for rank = 1:Rt
                    if tracksFrame(maskee).Rank <= Rt < tracksFrame(masker).Rank
                        masking(masker, maskee, :) = tracksFrame(masker).MagRes(:) - tracksFrame(maskee).MagRes(:);
                    else
                        masking(masker, maskee, :) = 0;
                    end
                end
                
            end
        end
    end
    
%%%%%%%%%%%%% MASKING SELECTION

    for masker = 1:size(tracks, 2)
        %Get max masking for the current masker
        Mask = selectMasking(masking(masker,:,:), nF); 
        
        %Smoothing
        for Fbin = 1:nF
            if(Mask(Fbin) ~= 0)
                tracksFrame(masker).filterM(Fbin) = EMA(Mask(Fbin), filter_d(Fbin).M, alpha);
                %For now we have 10 fixed filter centers
                %tracksFrame(masker).filterCenter = EMA(center(Fbin), filter_d(Fbin).center, alpha);
            else
                tracksFrame(masker).filterM(Fbin) = 0;
                %tracksFrame(masker).filterCenter(Fbin) = 0;
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%% FILTERING
    y = zeros(nF, size(tracksFrame(1).samples, 1));
    for track = 1:size(tracksFrame, 2)
        for i = 1:nF
            %[num,den] = iirpeak(filter_d(i).center/fs*2,(tracksFilt(i).center+tracksFilt(i).center/(Q*2))/fs*2);
            %[num,den] = iirpeak(0.0013,0.0014);
            %Y = filter(num, den, tracksFrame(track).samples);
            %y(i, :) =   Y; %tracksFrame(track).samples;%;%cv(1:M); %
        end
        %sorry about that
        [a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9), a(10)] = filter_d(1:10).center;
        
        %EQ   eq_filter(x, fc, Q, M, fs)
        y = eq_filter(tracksFrame(track).samples, a, zeros(1, 10)+Q, tracksFrame(track).filterM(1:10).*(-1), fs);
        
        %Sum the overlapping part of samples with new ones
        tracksFilt(track).samples(w_pos:w_pos+M-1) = [tracksFilt(track).samples(w_pos:w_pos+M*overlapR-1); zeros(M*(1-overlapR), 1)]' + y'; %[tracksFilt(track).samples, y];
    end

    % First frame alpha is 0.
    if(alpha == 0)
        alpha = exp(-1/(2*fs));
    end
    
    w_pos = w_pos + M*(1-overlapR);
    frame =  frame+1;
    fprintf("Frame %i/%i\n",frame, size(tracks(1).samples, 1)/M/(1-overlapR));
end
toc

%% Normalize
for i = 1:nTracks
    tracksFilt(i).samples = tracksFilt(i).samples/max(tracksFilt(i).samples)*0.8; %0.8 because want to be sure
end


%% WRITE AUDIOS
for i = 1:nTracks
    audiowrite("audios_rendered/FILT-"+tracks(i).name+".wav", real(tracksFilt(i).samples), fs);
end

for i = 1:nTracks %Rewrite original audios bc it could have been cut
    audiowrite("audios_rendered/ORIG-"+tracks(i).name+".wav", real(tracks(i).samples), fs);
end

%%
sound(real(tracksFilt(1).samples), fs)

%%
clear sound






