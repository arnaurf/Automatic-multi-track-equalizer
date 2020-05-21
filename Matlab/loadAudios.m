function audios = loadAudios(directory, winSize, overlapR, t0, tf)
    overlapSamples = ceil(overlapR*winSize);
    stepL = winSize - overlapSamples;
    files = ls(directory+"/*.wav"); %read files on the directory
    nTracks = size(files,1);
    audios = cell2table(cell(nTracks,3));
    audios.Properties.VariableNames = ["samples", "fs", "fileName"];

    if ~exist('t0','var') %If t0 doesn't exist, read from the begining
        t0 = 0;
    end
    
    j = 1;
    for i = 1:nTracks %reading all files on directory
        [x, fs] = audioread(directory+"/"+string(files(i,:)));

        if ~exist('tf', 'var') %if tf variable doesn't exist read the whole audio
            tf=size(x,1)/fs;
        end
        if size(x,1)/fs < tf || tf < 0 %check if tf is inside the audio length
            sprintf("End time is greater than audio length");
            return;
        end
        if t0 < 0 || t0 > size(x,1)/fs   %check if t0 is inside the audio length
            sprintf("Initial time must be greater than 0");
            return;
        end
        if tf < t0  
            sprintf("tf must be greater than t0")
            return;
        end

        x = x(floor(t0*fs+1):floor(tf*fs), :); %read only a section of the audio
        x(end: floor(size(x)/stepL)*stepL + winSize + 1, :) = 0; %Fill input with 0 to fit all windows
        
%         if(size(x, 2) > 1) %Stereo to mono
%             x = sum(x, 2)/size(x, 2);
%         end
        
        audios.samples(i,:) = {x'};
        audios.fs(i,:) = {fs};
        audios.fileName(i,:) = {files(i,files(i,:) ~= ' ')};
        j = j + 1;
    end    
end
