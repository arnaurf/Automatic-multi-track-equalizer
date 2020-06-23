/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <string>


//==============================================================================
MtequalizerAudioProcessorEditor::MtequalizerAudioProcessorEditor (MtequalizerAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{


	slider.setSliderStyle(Slider::LinearBarVertical);
	slider.setRange(0.0, 5.0, 0.001);
	slider.setTextBoxStyle(Slider::NoTextBox, false, 90, 0);
	slider.setPopupDisplayEnabled(true, false, this);
	slider.setTextValueSuffix("S");
	slider.setValue(1.0);
	slider.setSliderStyle(Slider::SliderStyle::Rotary);
	slider.setTextBoxStyle(Slider::TextBoxBelow, true, 70, 15);
	addAndMakeVisible(&slider);

	slider2.setSliderStyle(Slider::LinearBarVertical);
	slider2.setRange(0.0, 3.0, 0.01);
	slider2.setTextBoxStyle(Slider::NoTextBox, false, 90, 0);
	slider2.setPopupDisplayEnabled(true, false, this);
	slider2.setTextValueSuffix("Gain");
	slider2.setValue(1.0);
	slider2.setSliderStyle(Slider::SliderStyle::Rotary);
	slider2.setTextBoxStyle(Slider::TextBoxBelow, true, 70, 15);
	addAndMakeVisible(&slider2);

	test = false;

	addAndMakeVisible(&button1, 1);
	addAndMakeVisible(&plusBut, 1);
	addAndMakeVisible(&menusBut, 1);
	addAndMakeVisible(&label, 1);
	addAndMakeVisible(&label2, 1);

	button1.onClick = [this] { updateToggleState(&button1, "Male");   };
	button1.setRadioGroupId(1);
	//button1.setButtonText("Test");

	plusBut.onClick = [this] { updateToggleState(&plusBut, "+");   };
	plusBut.setRadioGroupId(2);
	plusBut.setButtonText("+");
	
	menusBut.onClick = [this] { updateToggleState(&menusBut, "-");   };
	menusBut.setRadioGroupId(2);
	menusBut.setButtonText("-");

	slider.addListener(this);
	slider2.addListener(this);

	button1.addListener(this);
	plusBut.addListener(this);
	menusBut.addListener(this);

	label.setFont(Font(14.0f, Font::plain));
	label.setJustificationType(juce::Justification::centred);

	label2.setFont(Font(14.0f, Font::plain));
	label2.setJustificationType(juce::Justification::centred);

	label.setText("nF ="+std::to_string(processor.aF), dontSendNotification);
	label2.setText("Toggle eq normalization", dontSendNotification);

	button1.setButtonText("off");

	// Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (400, 300);
}

MtequalizerAudioProcessorEditor::~MtequalizerAudioProcessorEditor()
{
}

void MtequalizerAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
	if(slider == &this->slider)
		processor.gain = slider->getValue();
	if (slider == &this->slider2)
		processor.post_gain = slider->getValue();

}
void MtequalizerAudioProcessorEditor::buttonClicked(Button* button) {
	if (button == &this->button1) {
		processor.eq_normalize = !processor.eq_normalize;
		String state = processor.eq_normalize ? "on" : "off";
		Logger::outputDebugString(" State changed to " + state);
		button->setButtonText(state);
	}
	else if (button == &this->plusBut) {
		if(processor.nF < 10)
			processor.nF += 1;
	}
	else if (button == &this->menusBut) {
		if (processor.nF > 3)
			processor.nF -= 1;
	}
	label.setText("nF =" + std::to_string(processor.nF), dontSendNotification);



}

void MtequalizerAudioProcessorEditor::updateToggleState(Button* button, String name) {
	/*
	auto state = button->getToggleState();
	String stateString = state ? "ON" : "OFF";
	button->setClickingTogglesState(true);

	Logger::outputDebugString(name + " Button changed to " + stateString);*/
}
//==============================================================================
void MtequalizerAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    g.setColour (Colours::white);
    //g.setFont (15.0f);
    //g.drawFittedText ("Hello World!", getLocalBounds(), Justification::, 1);

}

void MtequalizerAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	slider.setBounds(100, 10, 80, 80);

	slider2.setBounds(10, 10, 80, 80);


	// nF (+)(-)
	label.setBounds(185, 55, 100, 50);
	plusBut.setBounds(203, 37, 30, 30);
	menusBut.setBounds(237, 37, 30, 30);

	// Toggle maskee compensation
	label2.setBounds(280, 60, 100, 50);
	button1.setBounds(305, 20, 50, 50);

}
