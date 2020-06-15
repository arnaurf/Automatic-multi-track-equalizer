/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "Biquad.h"
#include "BiquadFilter.h"
#include <algorithm>
#include "utils.h"

#define MAX_CHANNELS 5


//==============================================================================
MtequalizerAudioProcessor::MtequalizerAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
         .withInput("Input1", AudioChannelSet::stereo(), true)
         .withInput("Input2", AudioChannelSet::stereo(), false)
         .withInput("Input3", AudioChannelSet::stereo(), false)
         .withInput("Input4", AudioChannelSet::stereo(), false)
         .withInput("Input5", AudioChannelSet::stereo(), false)
                      #endif
         .withOutput("Output1", AudioChannelSet::stereo(), true)
         .withOutput("Output2", AudioChannelSet::stereo(), false)
         .withOutput("Output3", AudioChannelSet::stereo(), false)
         .withOutput("Output4", AudioChannelSet::stereo(), false)
         .withOutput("Output5", AudioChannelSet::stereo(), false)
                     #endif
                       )
#endif
{


	last_gain = 1;
	gain = 1;
	eq_normalize = false;
	post_gain = 1;

    this->Q = 2;
    this->winSize = this->getBlockSize();
    this->overlapR = 0.25;
    this->overlapSamples = winSize * overlapR;
    this->stepL = winSize - overlapSamples;
    this->alpha = 0;
    this->Rt = 3;
    this->nF = 10;
    this->aF = 10;
    this->nTracks = this->getBusCount(true);
    this->firstFrame = true;


    filters = std::vector<std::vector<Filter>>(nTracks);
    for (int i = 0; i < nTracks; i++)
        filters[i] = std::vector<Filter>(nF);

    for (int i = 0; i < nTracks; i++) {
        for (int j = 0; j < nF; j++) {
            filters[i][j].center = 31.25 * pow(2, j);
            filters[i][j].gain = 0;
        }
    }


    M = std::vector<std::vector<std::vector<float>>>(nTracks);
    for (int i = 0; i < nTracks; i++) {
        M[i] = std::vector<std::vector<float>>(nTracks);
        for (int ii = 0; ii < nTracks; ii++) {
            M[i][ii] = std::vector<float>(aF);
        }
    }

   
}



MtequalizerAudioProcessor::~MtequalizerAudioProcessor()
{
}

//==============================================================================
const String MtequalizerAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool MtequalizerAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool MtequalizerAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool MtequalizerAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double MtequalizerAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int MtequalizerAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int MtequalizerAudioProcessor::getCurrentProgram()
{
    return 0;
}

void MtequalizerAudioProcessor::setCurrentProgram (int index)
{
}

const String MtequalizerAudioProcessor::getProgramName (int index)
{
    return {};
}

void MtequalizerAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void MtequalizerAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..



    this->alpha = 0;
    this->winSize = samplesPerBlock;
    firstFrame = true;

    // CLEAR ALL THE BUFFERS
    mPrevBuffers.clear();
    mPrevOverlapBuffers.clear();
    mCurrOverlapBuffers.clear();
    mCurrBuffers.clear();
    for (int i = 0; i < MAX_CHANNELS; i++) {
        mPrevBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mPrevBuffers[i].clear();

        mPrevOverlapBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mPrevOverlapBuffers[i].clear();

        mCurrOverlapBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mCurrOverlapBuffers[i].clear();

        mCurrBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mCurrBuffers[i].clear();
    }
    
    window = std::vector<double>(samplesPerBlock);
    hanning(window, samplesPerBlock);
    
}

void MtequalizerAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool MtequalizerAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else


    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
       && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif


//Returns on the MagRes array the magnitude of each frequency band inside of "buffer".
void MtequalizerAudioProcessor::getMagRes(float MagRes[], const AudioBuffer<float> buffer, bool isStereo) {

    //Init vars
    const float* channelData;
    AudioBuffer<float> aux = AudioBuffer<float>(buffer.getNumChannels(), buffer.getNumSamples());

    //If it is stereo, down to mono
	bool nonZero;
     for (int i = 0; i < buffer.getNumSamples(); i++) {
		 if (buffer.getSample(0, i) != 0)
			 nonZero = true;

		 if (isStereo) 
			aux.setSample(0, i, (buffer.getSample(0, i) + buffer.getSample(1, i))/2);
		 else
			aux.setSample(0, i, buffer.getSample(0, i));
	 }
		//channelData = buffer.getReadPointer(0);

	channelData = aux.getReadPointer(0);

    std::vector<double> xx;
    for (int i = 0; i < this->aF; i++) { //for each analisis frequency band

		if(!nonZero)
			MagRes[i] = -INFINITY;
		else {
			xx = filter2(&b[i * nTracks], &a[i * nTracks], channelData, 5, this->winSize);//bandpass filter it
			float rmsX = rms(xx, winSize); //calculate the RMS magnitude

			if (rmsX != 0)
				MagRes[i] = 20 * log10(rmsX);
			else
				MagRes[i] = -INFINITY;
		}
    }
    
}


//Returns inside the Rank array the ranking of the most important frequency bands of the signal, being 1 the most relevant and 10 the less relevant.
void MtequalizerAudioProcessor::getRank(float* Rank, float magnitude[]) {

    //Init vars
    std::vector<float> magnitude_cpy(magnitude, magnitude+aF);//
    float* max_value;


    for (int i = 0; i < this->aF; i++) {
        max_value = max(magnitude_cpy, this->aF); //Take the maximum value on magnitude_cpy. max_value = [value, index]
        magnitude_cpy[(int)max_value[1]] = -200;  //Set this value to -200 so it wont be the max value the next iteration
        Rank[(int)max_value[1]] = i;              //Set the rank array on that index to the iteration number (first max value found, is the most important)
		delete[] max_value;
    }
}

//x[5][10]. On a given track, it returns for each frequency band the bigger masking is producing to the other tracks.
void MtequalizerAudioProcessor::selectMasking(std::vector<float>& masking, std::vector<std::vector<float>> x) {

    std::vector<float> aux(aF);
    for (int iBand = 0; iBand < aF; iBand++) { //for each band (aF, analisisFilter bands)

        //Put all the masking of the other tracks on a single array.
        for (int iTrack = 0; iTrack < nTracks; iTrack++) {
            aux[iTrack] = x[iTrack][iBand];
        }
        float* max_value = max(aux, aF); //Then, the return the maximum value on that array.

        //Notice that as we computed the amount of masking as the difference between the magnitude of the masker and the maskee,
        //the masking amount is positive. But we want the "gain" compansation to be negative, that why (-1).
        masking[iBand] = max_value[0]*(-1);

		delete[] max_value;
    }
}

//if buffer does not contain any data, return false. If buffer contains 1 channel, return true.
//if buffer conaints 2 channels, but both contains the same data return true. Else, return false.
bool MtequalizerAudioProcessor::isMono(AudioBuffer<float> buffer) {

    if (buffer.getNumChannels() == 0)
        return false;

    if (buffer.getNumChannels() == 1)
        return true;

    //Check if we found differences between both channels. If it does, it is not mono.
    for (int i = 0; i < buffer.getNumSamples(); i++) {
        if (buffer.getSample(0, i) != buffer.getSample(1, i))
            return false;
    }
    return true;
}


//Equalize the track with the filter.
void MtequalizerAudioProcessor::eqfilter(Track& track, std::vector<Filter> filter, int Q) {

    Biquad lpFilter = Biquad();	// create a Biquad, lpFilter;
    for (int iBand = 1; iBand < this->nF; iBand++) {
       
        //If the gain is 0 (or too small) avoid filtering.
        if (filter[iBand].gain < -0.1) {

			
            BiquadFilter peakFilter = BiquadFilter(filter[iBand].center, Q, this->getSampleRate(), filter[iBand].gain*this->gain);


            //Equalize first channel
            int channel = 0;
            filter3(track.buffer->getWritePointer(channel), peakFilter, track.buffer->getReadPointer(channel), winSize);

            // CHECK FOR STEREO
            if (track.buffer->getNumChannels() == 2) { //If its stereo, equalize secondary channel
                channel = 1;
                //Reaper sometimes converts mono sources to stereo, so we have to check
                //If it is truly stereo, we equalize the second channel. If it a copy of the first channel, then we copy the previous result "Y"
                if (track.stereo) {
                    filter3(track.buffer->getWritePointer(channel), peakFilter, track.buffer->getReadPointer(channel), winSize);
                }
                else {
                    track.buffer->copyFrom(channel, 0, track.buffer->getReadPointer(0), winSize);
                }
            }

        }
        
    } 



     

}


//Thats the method that checks for the masked frequency bands of each track and equalize each track to reduce the masking.
void MtequalizerAudioProcessor::reduceMasking(std::vector<AudioBuffer<float>>& buffers) {
    
    if (buffers.size() <= 1 || nTracks <= 1) //If we only have one track enabled, nothing to compare with
        return;
    

    //Instert all track's info into a struct "Track" (MagRes, Ranking, Buffer and Stereo bool)
    std::vector<Track> track(nTracks);
    for (int iTrack = 0; iTrack < nTracks; ++iTrack) {

        track[iTrack].buffer = &buffers[iTrack];            //Save the buffer pointer
        track[iTrack].stereo = !isMono(buffers[iTrack]);    //Check if it is truly stereo
        getMagRes(track[iTrack].MagRes, buffers[iTrack], track[iTrack].stereo); // Save the magnitude of 10 frequency bands
        getRank(track[iTrack].Rank, track[iTrack].MagRes);  //Save the ranking of most important frequency bands

    }
    
    // Create Masking Matrix. It is, a comparation of which tracks are masking to others in a specific frequency band
    for (int masker = 0; masker < nTracks; ++masker) { //masker -> the track is creating the masking
        for (int maskee = 0; maskee < nTracks; ++maskee) { //maskee -> the track is being masked
            if (masker != maskee) { //don't compare with itself

                for (int i_band = 0; i_band < this->aF; i_band++) { //For each freq. band

                    if ((track[masker].Rank[i_band] < this->Rt) && (this->Rt <= track[maskee].Rank[i_band])
						&& track[masker].MagRes[i_band] > -50 && track[maskee].MagRes[i_band] > -50)
                        M[masker][maskee][i_band] = track[masker].MagRes[i_band] - track[maskee].MagRes[i_band];
                    else
                        M[masker][maskee][i_band] = 0;
                }

            }else { //if its comparing with itself, just put 0 masking.
                for (int i_band = 0; i_band < this->aF; i_band++)
                    M[masker][maskee][i_band] = 0;
            }
        }
    }

    
    //Now, for each track select the bigger amount of masking is producing to the others and smooth the change with the previous frame.
    std::vector< std::vector<float>> masking(nTracks);
	for (int i_masker = 0; i_masker < nTracks; i_masker++) {
		masking[i_masker] = std::vector<float>(nF);
		selectMasking(masking[i_masker], M[i_masker]);
	}
		
	if (eq_normalize) {
		for (int iB = 0; iB < nF; iB++) {
			int count = 0;
			for (int iT = 0; iT < nTracks; iT++) {
				if (masking[iT][iB] < 0)
					count++;
			}
		}
	}
	for (int i_masker = 0; i_masker < nTracks; i_masker++) {

        for (int i_bin = 0; i_bin < aF; i_bin++) {
            if (masking[i_masker][i_bin] != 0) {
                filters[i_masker][i_bin].gain = EMA(masking[i_masker][i_bin], filters[i_masker][i_bin].gain, 0.99f);
                //filters[i_masker][i_bin].gain = masking[i_bin];
            }
            else {
				filters[i_masker][i_bin].gain = EMA(0, filters[i_masker][i_bin].gain, 0.99f);
                //filters[i_masker][i_bin].gain = 0;
            }
        }

    }

    // Now we know the amount of masking, so we equalize each track.
    for (int iTrack = 0; iTrack < nTracks; iTrack++) {
        eqfilter(track[iTrack], filters[iTrack], Q);
    }


    //We use the alpha value for the smoothing process. Only the first frame alpha=0 (no smoothing)
    if (!alpha)
        alpha = exp(-1 / (2 * this->getSampleRate()));
 
    
}

//Thats is called each frame by the DAW.
void MtequalizerAudioProcessor::processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{

    auto totalNumInputChannels  = getTotalNumInputChannels();

    // Check which tracks are enabled and save them on mCurrBuffers
    this->nTracks = 0;
    int iTrack = 0;
    for (int i = 0; i < MAX_CHANNELS; i++) {
        Bus* bus = getBus(true, i);
        if (bus->isEnabled()) {
            this->nTracks += 1;
            mCurrBuffers[iTrack] = bus->getBusBuffer(buffer);
            iTrack += 1;
        }
        
    }
     this->nTracks = iTrack;  //Save the amount of actual enabled tracks
     
     
    //Complete mOverlapBuffers. This is the overlap part of the previous frame.
    for (int i = 0; i < nTracks; i++) {
        for (int channel = 0; channel < mCurrBuffers[i].getNumChannels(); channel++) {
            mCurrOverlapBuffers[i].copyFrom(channel, 0, mPrevBuffers[i], channel, winSize / 2, winSize / 2);
            mCurrOverlapBuffers[i].copyFrom(channel, winSize / 2, mCurrBuffers[i], channel, 0, winSize / 2);
        }
    }
    
    // Process the previous frame and the overlapping part of the current one. We process two "sub frames" each frame.
    reduceMasking(mPrevBuffers);
    reduceMasking(mCurrOverlapBuffers);

    
    // Construct the frame. We use the previous frame, plus the overlapping part of the pre-previous frame and the overlapping part of the current frame
    for (int i = 0; i < nTracks; i++) {
        for (int channel = 0; channel < mCurrBuffers[i].getNumChannels(); channel++) {

            float aux;
            for (int s = 0; s < winSize; s++) { //Iterate all samples
                aux = mCurrBuffers[i].getSample(channel, s); //save it for later

                //The first half of the current frame is the previous frame + the pre-previous frame overlapping part
                if (s < winSize / 2)  
                    mCurrBuffers[i].setSample(channel, s, mPrevOverlapBuffers[i].getSample(channel, s + winSize / 2) * window[int(s + winSize / 2)] + mPrevBuffers[i].getSample(channel, s) * window[s]);

                //The second half of the current frame, is the previous frame + the current frame overlapping part
                else   
                    mCurrBuffers[i].setSample(channel, s, mCurrOverlapBuffers[i].getSample(channel, s - winSize / 2) * window[int(s - winSize / 2)] + mPrevBuffers[i].getSample(channel, s) * window[s]);

                // For the next frame, the previous frame is the current one now.
                mPrevBuffers[i].setSample(channel, s, aux);
            }
            
            //The pre-previous overlapping part on the next frame is the current overlapping buffer.
            mPrevOverlapBuffers[i].copyFrom(channel, 0, mCurrOverlapBuffers[i], channel, 0, winSize);

        }

    }
    

	for (int i = 0; i < nTracks; i++) {
		if (this->last_gain != this->post_gain) {
			mCurrBuffers[i].applyGainRamp(0, this->winSize, last_gain, post_gain);
			last_gain = post_gain;
		}else{
			mCurrBuffers[i].applyGain(post_gain);
			//apply gain
		}

	}
    
    //Clear midi messages, we don't produce any.
    midiMessages.clear();
}

//==============================================================================
bool MtequalizerAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* MtequalizerAudioProcessor::createEditor()
{
    return new MtequalizerAudioProcessorEditor (*this);
}

//==============================================================================
void MtequalizerAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void MtequalizerAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new MtequalizerAudioProcessor();
}
