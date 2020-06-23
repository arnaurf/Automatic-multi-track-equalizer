/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class MtequalizerAudioProcessorEditor : public AudioProcessorEditor, private Slider::Listener, private Button::Listener
{
public:
    MtequalizerAudioProcessorEditor (MtequalizerAudioProcessor&);
    ~MtequalizerAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;

private:

	enum RadioButtonIds
	{
		GenderButtons = 1001
	};

    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
	void sliderValueChanged(Slider* slider) override;
	void updateToggleState(Button* button, String name);
	void buttonClicked(Button* button) override;


    MtequalizerAudioProcessor& processor;
	

	TextButton button1, plusBut, menusBut;
	Label label, label2;
	bool test;
	Slider slider, slider2, slider3;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtequalizerAudioProcessorEditor)
};
