/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once
#define B_COEF { 1.61102657579673e-05,0,-3.22205315159346e-05,0,1.61102657579673e-05, 3.73252044940333e-05, 0, -7.46504089880666e-05, 0, 3.73252044940333e-05,0.000148021987085363, 0, -0.000296043974170726, 0, 0.000148021987085363,0.000582076134415265, 0, -0.00116415226883053, 0, 0.000582076134415265,0.00225158265688761, 0, -0.00450316531377522, 0, 0.00225158265688761,0.00844269292908856, 0, -0.0168853858581771, 0, 0.00844269292908856,0.02995458220809, 0, -0.05990916441618, 0, 0.02995458220809,0.097631072937818, 0, -0.195262145875636, 0, 0.097631072937818,0.292893218813452, 0, -0.585786437626904, 0, 0.292893218813452,0.465036673649006, -0, -0.930073347298012, -0, 0.465036673649006 }
#define A_COEF { 1,-3.98861309864167,5.96590588446702,-3.96597246088351,0.988679675059301,1, -3.98251213329087, 5.94782014192403, -3.94810272349167, 0.982794719299833,1, -3.96476254351518, 5.89541975179491, -3.89654259839845, 0.965885460568819,1, -3.92850150864099, 5.79001068662419, -3.79444280174679, 0.932934731756612,1, -3.85308702992494, 5.57711422446639, -3.59437753745536, 0.870367477456467,1, -3.69181603431387, 5.14559727297041, -3.2110717379791, 0.757546944478831,1, -3.33500858117699, 4.27709277949262, -2.512537824396, 0.574061915083955,1, -2.52704741345686, 2.62123443505209, -1.38208812331399, 0.333333333333333,1, -0.732050807568877, 0.156961002899234, -0.125600061886466, 0.171572875253809,1, 1.37942528737802, 6.12074755323988e-08, -0.139102623328703, 0.240322786464271 }

#include <JuceHeader.h>

//==============================================================================
/**
*/
class MtequalizerAudioProcessor  : public AudioProcessor
{
public:

    struct sfilter {
        float center;
        float gain;
    }typedef Filter;

    struct strack {
        float MagRes[10];
        float Rank[10];
        AudioBuffer<float>* buffer;
        bool stereo;
    }typedef Track;

    float winSize, overlapR, overlapSamples, stepL, alpha = 0;
    int Rt, nF, aF, nTracks = 0;
    int Q = 0;
    double b[50] = B_COEF;
    double a[50] = A_COEF;
    std::vector<std::vector<Filter>> filter;
    std::vector<std::vector<std::vector<float>>> M;
    bool firstFrame;
    std::vector<double> window;

    //AudioBuffer<float>* buffers;
    std::vector<AudioBuffer<float>> mCurrBuffers; //current buffers, real i/o bufers
    std::vector<AudioBuffer<float>> mPrevBuffers; //last buffer (raw, unprecessed)
    std::vector<AudioBuffer<float>> mCurrOverlapBuffers;  //second half of the last frame buffer, and first half of the current buffer
    std::vector<AudioBuffer<float>> mPrevOverlapBuffers;

    void getMagRes(float MagRes[], const AudioBuffer<float> buffer, bool isStereo);
    void getRank(float* Rank, float magnitude[]);
    void selectMasking(std::vector<float> &masking, std::vector<std::vector<float>> x);
    void eqfilter(Track &track, std::vector<Filter> filter, int Q);
    void reduceMasking(std::vector<AudioBuffer<float>>& buffers);

    //==============================================================================
    MtequalizerAudioProcessor();
    ~MtequalizerAudioProcessor();

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    bool isMono(AudioBuffer<float> buffer);

private:
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtequalizerAudioProcessor)
};
