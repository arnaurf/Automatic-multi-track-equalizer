clear all; 
%% %%%%%%% INIT VARIABLES
% WINDOW & OVERLAP
winSize = 1024; %window size
overlapR = 0; %overlap %
overlapSamples = ceil(overlapR*winSize);
stepL = winSize - overlapSamples;

% LOAD AUDIOS -> loadAudios(directory, winSize, overlapR, t0, tf)
x_original = loadAudios("audios", winSize, overlapR); %t0 and tf optionals

nTracks = size(x_original, 1);
Length = size(x_original.samples{1},2);
fs = x_original.fs{1};


% for i=1:nTracks
%     t = x_original.samples{i};
%     x_original.samples{i} = x_original.samples{i}./max(x_original.samples{i}, [], 'all');
% end

y_filtered = cell(nTracks,1);
for i=1:nTracks
    y_filtered{i} = zeros(size(x_original.samples{i})); %create empty output
end
    
    
% INIT FILTER
Q = 2;
nF = 5; %number of filters EQ
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
W = hanning(winSize);
alpha = 0; %alpha value for EMA smoothing
S = 1;
eq_normalize = false;
weighted_mean = false;
variable_nF = false;

% EMPTY FRAME VARIABLE
frame = table(zeros(nTracks, winSize), zeros(nTracks, aF), zeros(nTracks, aF));
frame.Properties.VariableNames = ["Samples", "MagRes", "Rank"];

frameSamples = cell(nTracks,1);
frameMagnitudes = zeros(ceil(Length/stepL), nTracks, aF);
frameOverallMagnitudes = zeros(ceil(Length/stepL), nTracks);
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
       
        %%%%%%%%%%%% FEATURES EXTRACTION
        frameMagnitudes(nframe,iTrack,:) = getMagRes(samples_mono, b, a);        %get Magnitude of each band
        frameOverallMagnitudes(nframe, iTrack) = rms(samples_mono);
    end
    
    winPos = winPos + stepL;
    nframe =  nframe+1;
end   
%%
Magnitudes = zeros(nTracks, aF);
Rank = zeros(nTracks, aF);




for iT = 1:nTracks
    if weighted_mean == true
        weights = frameOverallMagnitudes(:, iT)./sum(frameOverallMagnitudes(:, iT)); % Compute weighted mean, depending on the magnitude of each frame
        mean_mag = sum(squeeze(frameMagnitudes(:, iT, :)) .* weights);
    else
        mean_mag = mean(squeeze(frameMagnitudes(:, iT, :)),1);
    end
    
    Magnitudes(iT,:) = mean_mag;%mean(squeeze(frameMagnitudes(:, iT, :)),1);
    Rank(iT,:) = getRank(Magnitudes(iT,:), aF)';   
end

%%%%%%%%%%%%%%%%%%%%% MASKING DETECTION
M = cell(nTracks);
for i_masker= 1:nTracks %for each masker
    for i_maskee = 1:nTracks %compare it with each possible maskee
        if i_masker ~= i_maskee %avoid comparing with itself
            for i_band = 1:aF %per each band
                if (Rank(i_maskee,i_band) <= Rt) && (Rt < Rank(i_masker, i_band)) %Equation (1) on reference [1]
                    M{i_masker, i_maskee}(i_band) = Magnitudes(i_masker,i_band) - Magnitudes(i_maskee,i_band);
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
Mask = zeros(nTracks, aF);
for i_masker = 1:nTracks %for each masker
    %Get max masking for the current masker on each band
    Mask(i_masker,:) = selectMasking(cell2mat(M(i_masker,:)'), nF); 
end

if eq_normalize == true
    for iBand = 1:aF
        Mask(:,iBand) = Mask(:,iBand)./sum(Mask(:,iBand)>0);
    end
end
%%%%%%%%%%%%%%%% FILTERING
for i_track = 1:nTracks
    %EQUALIZE   eq_filter(x, fc, Q, winSize, fs)
    y_filtered{i_track} = eq_filter(x_original.samples{i_track}, filter.center(i_track,:), zeros(1, 10)+Q, Mask(i_track,:).*(-1)*S, fs);
end

t = toc;
fprintf("Elapsed time is %f seconds, performance: %fx\n", t, Length/fs/t);
%% Normalize
% for i = 1:nTracks
%     tracksFilt(i).samples = tracksFilt(i).samples/max(tracksFilt(i).samples)*0.8; %0.8 because want to be sure
% end


%% WRITE AUDIOS
y_mixdown = 0;
for i = 1:nTracks
    audiowrite("audios_rendered/FILT-"+x_original.fileName{i}, real(y_filtered{i})', fs);
    y_mixdown = y_mixdown + y_filtered{i};
end


x_mixdown = 0;
for i = 1:nTracks %Rewrite original audios bc it could have been cut
    audiowrite("audios_rendered/ORIG-"+x_original.fileName{i}, real(x_original.samples{i,:}'), fs);
   x_mixdown = x_mixdown + x_original.samples{i};
end

% gain_compensation = rms(x_mixdown(1,:))/rms(y_mixdown(1,:));
% y_mixdown = y_mixdown.*gain_compensation;

x_mixdown = x_mixdown./max(x_mixdown, [], 'all').*0.99;
y_mixdown = y_mixdown./max(y_mixdown, [], 'all').*0.99;

audiowrite("audios_rendered/mixdownORIG-"+x_original.fileName{i}, real(x_mixdown)', fs);
audiowrite("audios_rendered/mixdownFILT-"+x_original.fileName{i}, real(y_mixdown)', fs);


%%
sound(real(y_filtered{1}), fs);

%%
clear sound


%%
smr_x = mean(smr(x_mixdown', fs));
smr_y = mean(smr(y_mixdown', fs));



% masking_table = zeros(1, 2);
% for maskee = 1:nTracks
%     
%     aux = 0;
%     for masker = 1:nTracks
%         if masker ~= maskee
%             masking = integratedLoudness(y_filtered{maskee}'+y_filtered{masker}',fs)/integratedLoudness(y_filtered{maskee}',fs)*100;
%         end
%     end
%     integratedLoudness(y_filtered{1}'+y_filtered{2}',fs)/integratedLoudness(y_filtered{1}',fs)*100;
%     
%     
% end

%% %%%%%%%%%%%%  REFERENCES  %%%%%%%%%%%%%
% [1] Hafezi, S., & Reiss, J. (2015). Autonomous Multitrack Equalization Based on Masking Reduction. Journal of the Audio Engineering Society, 63(5), 312–323. doi: 10.17743/jaes.2015.0021
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
