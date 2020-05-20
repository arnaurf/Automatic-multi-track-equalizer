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


int offset(int x, int y, int z, int xSize, int ySize) {
    return (z * xSize * ySize) + (y * xSize) + x;
}

float EMA(float x, float y0, float alpha) {
    return (1 - alpha) * x + alpha * y0;
}

void hanning(double* buffer, int size) {
    for (int i = 0; i < size; i++) {
        buffer[i] =  0.5 * (1 - cos(2 * M_PI * i / size));
    }
}

double* filter2(double b[], double a[], const float* X, int sizeF, int sizeX) {

    //float* z = (float*)calloc(sizeX, sizeof(float));
    for (int i = 0; i < sizeF; i++) {
        b[i] = b[i] / a[0];
        a[i] = a[i] / a[0];
    }

    double* Y = (double*)calloc(sizeX, sizeof(double));
    
    for (int n = 0; n < sizeX; n++) {
        double auxX = 0;
        double auxY = 0;
        for (int m = 0; m < sizeF; m++) {
            if(n-m >= 0)
                auxX = auxX + X[n - m] * b[m];
        }
        for (int m = 1; m < sizeF; m++) {
            if (n - m>= 0)
                auxY = auxY + Y[n - m] * a[m];
        }
        Y[n] = (auxX - auxY) / double(a[0]);
    }
    
    /*for (int m = 1; m < sizeX; m++) {
        Y[m] = b[0] * X[m] + z[0];
        for (int i = 2; i < sizeF; i++) {
            z[i - 1] = b[i] * X[m] + z[i] - a[i] * Y[m];
        }
    }
    //z = z[1:n - 1];*/
    return Y;
}

double* filter3(double b[], double a[], const float* X, int sizeB, int sizeA, int sizeX) {

    //float* z = (float*)calloc(sizeX, sizeof(float));
    for (int i = 0; i < sizeB; i++) {
        b[i] = b[i] / a[0];
    }
    for (int i = 0; i < sizeA; i++) {
        a[i] = a[i] / a[0];
    }
    double* Y = (double*)calloc(sizeX, sizeof(double));

    for (int n = 0; n < sizeX; n++) {
        double auxX = 0;
        double auxY = 0;
        for (int m = 0; m < sizeB; m++) {
            if (n - m >= 0)
                auxX = auxX + X[n - m] * b[m];
        }
        for (int m = 1; m < sizeA; m++) {
            if (n - m >= 0)
                auxY = auxY + Y[n - m] * a[m];
        }
        Y[n] = (auxX - auxY) / double(a[0]);
    }

    return Y;
}

float rms(double x[], int n)
{
    double sum = 0;

    for (int i = 0; i < n; i++)
        sum += pow(x[i], 2);

    return sqrt(sum / n);
}

float* max(float x[], int size) {

    float temp_max[2] = { x[0],0 };
    for (int i = 0; i < size; i++) {
        if (x[i] > temp_max[0]) {
            temp_max[0] = x[i];
            temp_max[1] = i;
        }
    }
    return temp_max;
}



//==============================================================================
MtequalizerAudioProcessor::MtequalizerAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
         .withInput("Input1", AudioChannelSet::stereo(), true)
         .withInput("Input2", AudioChannelSet::stereo(), true)
         .withInput("Input3", AudioChannelSet::stereo(), true)
         .withInput("Input4", AudioChannelSet::stereo(), true)
         .withInput("Input5", AudioChannelSet::stereo(), true)
                      #endif
         .withOutput("Output1", AudioChannelSet::stereo(), true)
         .withOutput("Output2", AudioChannelSet::stereo(), true)
         .withOutput("Output3", AudioChannelSet::stereo(), true)
         .withOutput("Output4", AudioChannelSet::stereo(), true)
         .withOutput("Output5", AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{

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


    //filter = (Filter*)malloc(sizeof(Filter) * nTracks * nF);
    filter = new Filter*[nTracks];
    for (int i = 0; i < nTracks; i++)
        filter[i] = new Filter[nF];


    //buffers = new AudioBuffer<float>[nTracks];
    

    
    float x[400];
    for (int i = 0; i < 400; i++) {
        x[i] = i * 0.01;
    }
    //float x[10] = { 1,2,3,4,5,6,7,8,9,10 };
    //double b[] = { 1 };
    //double a[] = { 1, 0.2 };
    double* y;

   //y = filter2(b, a, x, 5, 400);


    for (int i = 0; i < nTracks; i++) {
        for (int j = 0; j < nF; j++) {
            filter[i][j].center = 31.25 * pow(2, j);
            filter[i][j].gain = 0;
        }
    }
   // M = (float***)malloc(sizeof(float**) * nTracks);
    M = new float**[nTracks];
    for (int i = 0; i < nTracks; i++) {
        M[i] = new float* [nTracks];
        for (int ii = 0; ii < nTracks; ii++) {
            M[i][ii] = new float[aF];
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
   
    firstFrame = true;
    mPrevBuffers.clear();
    mPrevOverlapBuffers.clear();
    mCurrOverlapBuffers.clear();
    mCurrBuffers.clear();

    for (int i = 0; i < nTracks; i++) {
        mPrevBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mPrevBuffers[i].clear();

        mPrevOverlapBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mPrevOverlapBuffers[i].clear();

        mCurrOverlapBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mCurrOverlapBuffers[i].clear();

        mCurrBuffers.push_back(AudioBuffer<float>(2, samplesPerBlock));
        mCurrBuffers[i].clear();
    }
    window = new double[winSize];
    hanning(window, winSize);
    
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
void MtequalizerAudioProcessor::getMagRes(float MagRes[], const float* channelData) {
    int sizeF = sizeof(b) / sizeof(*b) / nTracks;
    double* xx = (double*)malloc(sizeof(double) * sizeF);
    for (int i = 0; i < sizeF; i++) {
        xx = filter2(&b[i * nTracks], &a[i * nTracks], channelData, 5, this->winSize);
        float rmsX = rms(xx, winSize);
        if (rmsX != 0)
            MagRes[i] = 20 * log10(rmsX);
        else
            MagRes[i] = -INFINITY;
    }

}

void MtequalizerAudioProcessor::getRank(float* Rank, float magnitude[]) {
    float* magnitude_cpy = (float*)malloc(sizeof(float) * this->aF);

    memcpy(magnitude_cpy, magnitude, sizeof(float) * this->aF);
    float* max_value;
    for (int i = 0; i < this->aF; i++) {
        max_value = max(magnitude_cpy, this->aF);
        magnitude_cpy[(int)max_value[1]] = -200;
        Rank[(int)max_value[1]] = i;
    }
}

//x[5][10]
void MtequalizerAudioProcessor::selectMasking(float* masking, float** x) {

    float* aux = new float[aF];
    for (int i = 0; i < aF; i++) { //for each band (aF)

        for (int j = 0; j < nTracks; j++) {
            aux[j] = x[j][i];
        }
        float* output = max(aux, aF);
        masking[i] = output[0]*(-1);
    }
}

void MtequalizerAudioProcessor::eqfilter(AudioBuffer<float>* buffer, float*input, Filter* filter, int Q) {

    Biquad* lpFilter = new Biquad();	// create a Biquad, lpFilter;
    BiquadFilter aux3 = BiquadFilter(100, 2, 48000, 2);
    for (int iBand = 0; iBand < this->nF; iBand++) {
        lpFilter->setBiquad(bq_type_peak, filter[iBand].center / this->getSampleRate(), Q, filter[iBand].gain);
        double* aux2 = lpFilter->getCoef();
        for (int channel = 0; channel < 2; channel++) {
            float * Y = (float*)filter3(&aux2[3], aux2, input, 3,2, winSize);
            *(buffer->getWritePointer(channel)) = *Y;

            /*for (int idx = 0; idx < this->winSize; idx++) {
                float aux = lpFilter->process(input[idx]);
                float a = 0;
                (buffer->getWritePointer(channel))[idx] = 1.0f * aux;
            }*/


        }
    } 



     

}

void MtequalizerAudioProcessor::reduceMasking(std::vector<AudioBuffer<float>>& buffers) {
    Track* track;
    track = new Track[nTracks];

    for (int channel = 0; channel < nTracks; ++channel) {
        track[channel].Samples = buffers[channel].getWritePointer(0);
        getMagRes(track[channel].MagRes, buffers[channel].getReadPointer(0));
        getRank(track[channel].Rank, track[channel].MagRes);
    }//IT WORKS!

    for (int masker = 0; masker < nTracks; ++masker) {
        for (int maskee = 0; maskee < nTracks; ++maskee) {
            if (masker != maskee) {
                for (int i_band = 0; i_band < this->aF; i_band++) {
                    if ((track[masker].Rank[i_band] < this->Rt) && (this->Rt <= track[maskee].Rank[i_band]))
                        M[masker][maskee][i_band] = track[masker].MagRes[i_band] - track[maskee].MagRes[i_band];
                    else
                        M[masker][maskee][i_band] = 0;
                }
            }
            else {
                for (int i_band = 0; i_band < this->aF; i_band++)
                    M[masker][maskee][i_band] = 0;
            }
        }
    }

    float* masking = (float*)malloc(sizeof(float) * aF);
    for (int i_masker = 0; i_masker < nTracks; i_masker++) {
        selectMasking(masking, M[i_masker]);

        for (int i_bin = 0; i_bin < aF; i_bin++) {
            if (masking[i_bin] != 0) {
                //filter[i_masker][i_bin].gain = EMA(masking[i_bin], filter[i_masker][i_bin].gain, this->alpha);
                filter[i_masker][i_bin].gain = masking[i_bin];
            }
            else {
                filter[i_masker][i_bin].gain = 0;
            }
        }

    }

    for (int iTrack = 0; iTrack < nTracks; iTrack++) {
        eqfilter(&buffers[iTrack], track[iTrack].Samples, filter[iTrack], Q);

        for (int idx = 0; idx < this->winSize; idx++) {
            for (int channel = 0; channel < 2; channel++) {
                (buffers[iTrack].getWritePointer(channel))[idx] = track[iTrack].Samples[idx];
                //(buffers[iTrack]->getWritePointer(channel))[idx] = track[iTrack].Samples[idx];
            }
        }
    }

    if (!alpha)
        alpha = exp(-1 / (2 * this->getSampleRate()));
}


void MtequalizerAudioProcessor::processBlock(AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    //auto totalNumOutputChannels = getTotalNumOutputChannels();   
    int totalChannels = 2;

    this->winSize = buffer.getNumSamples();


    //float* M = (float*)malloc(sizeof(float) * nTracks * nTracks * this->aF);

    //AudioBuffer<float>* buffers = (AudioBuffer<float>*) calloc(sizeof(AudioBuffer<float>), nTracks);
     
    AudioBuffer<float> aux;
     for (int i = 0; i < nTracks; i++) {
        Bus* bus = getBus(true, i);
        mCurrBuffers[i] = (bus->getBusBuffer(buffer));
    }
    
 
    //Complete mOverlapBuffers
    for (int i = 0; i < nTracks; i++) {
        for (int channel = 0; channel < totalChannels; channel++) {
            mCurrOverlapBuffers[i].copyFrom(channel, 0, mPrevBuffers[i], channel, winSize / 2, winSize / 2);
            mCurrOverlapBuffers[i].copyFrom(channel, winSize / 2, mCurrBuffers[i], channel, 0, winSize / 2);
        }
    }
   
    reduceMasking(mCurrBuffers);
    reduceMasking(mCurrOverlapBuffers);

    
    for (int i = 0; i < nTracks; i++) {

        for (int channel = 0; channel < totalChannels; channel++) {

            float aux;
            /*for (int s = 0; s < winSize; s++) {
                
                mCurrBuffers[i].setSample(channel, s, mPrevBuffers[i].getSample(channel, s));
                mPrevBuffers[i].setSample(channel, s, aux);
            }
            */
            for (int s = 0; s < winSize; s++) { //Iterate all samples

                aux = mCurrBuffers[i].getSample(channel, s);
                if (s < winSize / 2)  //The first half, output = LastBuffer
                    mCurrBuffers[i].setSample(channel, s, mPrevOverlapBuffers[i].getSample(channel, s + winSize / 2) * window[int(s + winSize / 2)] + mPrevBuffers[i].getSample(channel, s) * window[s]);

                else   //The second half, output = second half of LastBuffer + first half of mOverlapBuffer
                    mCurrBuffers[i].setSample(channel, s, mCurrOverlapBuffers[i].getSample(channel, s - winSize / 2) * window[int(s - winSize / 2)] + mPrevBuffers[i].getSample(channel, s) * window[s]);

                mPrevBuffers[i].setSample(channel, s, aux);
            }
            

            mPrevOverlapBuffers[i].copyFrom(channel, 0, mCurrOverlapBuffers[i], channel, 0, winSize);

        }

    }

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
