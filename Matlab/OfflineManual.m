clear all; 
%% %%%%%%% INIT VARIABLES

% ============ USER PARAMETERS =========== %
eq_normalize = false;
weighted_mean = true;
S = 2; %Scale for the filter gain  (gain*S)
Rt = 3; %max Rank
Q = 2;
nF = 5; %number of filters EQ
% ======================================== %





% =============== PROGRAM ================ %
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


y_filtered = cell(nTracks,1);
for i=1:nTracks
    y_filtered{i} = zeros(size(x_original.samples{i})); %create empty output
end
    
% INIT FILTERS

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

W = hanning(winSize);

%% =====MANUAL ESSENTIAL-NONESSENTIAL=====
band_classif = zeros(nTracks, aF);
for iT=1:nTracks
    figure(1);
    showSpectrum(x_original.samples{iT}, fs, fc, string(x_original.fileName{iT}));
    loop = 1;
    
    while loop
        inp = input("\n - Enter " + string(Rt) + " most essential bands for track \'"+string(x_original.fileName{iT})+"\' separated by ',' (example: \'1,2,9\') (see bands and spectrum in the figure)\n", 's');
        inp_num = str2double(inp(inp~=','));
        inp_char = inp(inp~=',');
        if ~isnan(inp_num) && length(inp_char) == Rt && length(unique(inp_char)) == Rt
            loop = 0;
            for i=1:Rt
                band_classif(iT, str2double(inp_char(i))+1) = 1;
            end
        else 
            fprintf("Invalid input, try again"); 
        end
    end
    
    
end
%%========================================


%% EMPTY FRAME VARIABLE
frame = table(zeros(nTracks, winSize), zeros(nTracks, aF), zeros(nTracks, aF));
frame.Properties.VariableNames = ["Samples", "MagRes", "Rank"];

frameSamples = cell(nTracks,1);
frameMagnitudes = zeros(floor(Length/stepL), nTracks, aF);
frameOverallMagnitudes = zeros(floor(Length/stepL), nTracks);

% MAIN LOOP: for each frame of size winSize: Compute magnitudes an rankings for each frame
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
        if any(samples_mono) %If the frame is not empty
            frameMagnitudes(nframe,iTrack,:) = getMagRes(samples_mono, b, a);        %get Magnitude of each band
            frameOverallMagnitudes(nframe, iTrack) = rms(samples_mono);
        else
            frameMagnitudes(nframe,iTrack,:) = -Inf;
            frameOverallMagnitudes(nframe, iTrack) = 0;
        end
        
        
    end
    
    winPos = winPos + stepL;
    nframe =  nframe+1;
end   

%%
% ====== Compute masking and filters for each track ====== %

Magnitudes = zeros(nTracks, aF);
Rank = zeros(nTracks, aF);
for iT = 1:nTracks
    %Get magnitudes for that track
    trackMagnitudes = squeeze(frameMagnitudes(:, iT, :));
    trackMagnitudes = 10.^(trackMagnitudes./20); %Transform to linear for the operations
    
    %Take into account the rms of each frame for the magnitude mean 
    if weighted_mean == true
        weights = frameOverallMagnitudes(:, iT)./sum(frameOverallMagnitudes(:, iT)); % Compute weighted mean, depending on the magnitude of each frame
        mean_mag = sum(trackMagnitudes.* weights);
    else
        mean_mag = mean(trackMagnitudes,1); %Normal mean
    end
    
    Magnitudes(iT,:) = 20*log10(mean_mag); %Transform to logarithmic
    Rank(iT,:) = getRank(Magnitudes(iT,:), aF)';   
end

% ====================== MASKING DETECTION ====================== %
M = cell(nTracks);
for i_masker= 1:nTracks %for each masker
    for i_maskee = 1:nTracks %compare it with each possible maskee
        if i_masker ~= i_maskee %avoid comparing with itself
            for i_band = 1:aF %per each band
                if band_classif(i_maskee,i_band)==1 && band_classif(i_masker, i_band)==0 ...
                        && Magnitudes(i_masker,i_band) ~= -Inf && Magnitudes(i_maskee,i_band) ~= -Inf %Equation (1) on reference [1]
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
% ===========================================================
    
% ==================== MASKING SELECTION ====================
Mask = zeros(nTracks, aF);
for i_masker = 1:nTracks %for each masker
    %Get max masking for the current masker on each band
    Mask(i_masker,:) = selectMasking(cell2mat(M(i_masker,:)')*(-1), nF); % INPUT MUST BE NEGATIVE
end

% Allowing: Don't allow to apply too much filtering on the same band at overall
if eq_normalize == true
    for iBand = 1:aF
        Mask(:,iBand) = Mask(:,iBand)./sum(Mask(:,iBand)>0);
    end
end

% ====================== FILTERING ====================== %
for i_track = 1:nTracks
    %EQUALIZE   eq_filter(x, fc, Q, winSize, fs)
    y_filtered{i_track} = eq_filter(x_original.samples{i_track}, filter.center(i_track,:), zeros(1, 10)+Q, Mask(i_track,:).*(-1)*S, fs);
end
% ======================================================= %

t = toc;
fprintf("Elapsed time is %f seconds, performance: %fx\n", t, Length/fs/t);


%% ====================== WRITE AUDIOS ====================== %%
% Filtered signal
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

audiowrite("audios_rendered/mixdownORIG-"+x_original.fileName{i}, real(x_mixdown)', fs);
audiowrite("audios_rendered/mixdownFILT-"+x_original.fileName{i}, real(y_mixdown)', fs);


%% ============================ EVALUATION ================================
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

%% EVALUATION
iTrack_compare = 6;
T1 = mask_curves{iTrack_compare,1};
T2 = mask_curves{iTrack_compare,2};
ind1 = zeros(1, size(T1,2));
ind2 = zeros(1, size(T2,2));
[ind1(1,:), ind1(2,:)] = ismember(T1(2,:),T2(2,:));
[ind2(1,:), ind2(2,:)] = ismember(T2(2,:),T1(2,:));
plot(T2(2,ind1(2,ind1(2,:)~=0)), T2(1,ind1(2,ind1(2,:)~=0)) - T1(1,ind2(2,ind2(2,:)~=0)) );
grid on;
%% %%%%%%%%%%%%  REFERENCES  %%%%%%%%%%%%%
% [1] Hafezi, S., & Reiss, J. (2015). Autonomous Multitrack Equalization Based on Masking Reduction. Journal of the Audio Engineering Society, 63(5), 312–323. doi: 10.17743/jaes.2015.0021
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
