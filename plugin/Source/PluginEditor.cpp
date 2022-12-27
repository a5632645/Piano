#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PianoAudioProcessorEditor::PianoAudioProcessorEditor (PianoAudioProcessor& p_)
    : ProcessorEditor (p_), proc (p_)
{
    additionalProgramming = "Clayton Otey";
    setName ("main");
    
    keyboard.setName ("keys");
    keyboard.setOpaque (false);
    keyboard.setMidiChannel (1);
    keyboard.setMidiChannelsToDisplay (1 << 0);
    keyboard.setKeyWidth (20);
    keyboard.setAvailableRange (PIANO_MIN_NOTE, PIANO_MAX_NOTE);
    keyboard.setScrollButtonsVisible (false);
    addAndMakeVisible (keyboard);

    setGridSize (15, 5);
}

PianoAudioProcessorEditor::~PianoAudioProcessorEditor()
{
}

//==============================================================================
void PianoAudioProcessorEditor::paint (juce::Graphics& g)
{
    ProcessorEditor::paint (g);
}

void PianoAudioProcessorEditor::resized()
{
    ProcessorEditor::resized ();
}
