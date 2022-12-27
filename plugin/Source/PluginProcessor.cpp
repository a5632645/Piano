#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PianoAudioProcessor::PianoAudioProcessor()
{
    midiOut.ensureSize (1024);
}

PianoAudioProcessor::~PianoAudioProcessor()
{
}

//==============================================================================
void PianoAudioProcessor::stateUpdated()
{
}

void PianoAudioProcessor::updateState()
{
}

//==============================================================================
void PianoAudioProcessor::reset()
{
    Processor::reset();
}

void PianoAudioProcessor::prepareToPlay (double newSampleRate, int newSamplesPerBlock)
{
    Processor::prepareToPlay (newSampleRate, newSamplesPerBlock);

    piano.init (float (newSampleRate), newSamplesPerBlock);
}

void PianoAudioProcessor::releaseResources()
{
}

void PianoAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midi)
{
    juce::ScopedNoDenormals noDenormals;

	buffer.clear ();
    auto numSamples = buffer.getNumSamples();

	keyState.processNextMidiBuffer (midiOut, 0, numSamples, true);

    auto ptr = (float**)buffer.getArrayOfWritePointers();
    piano.process (ptr, ptr, buffer.getNumSamples(), midi);
}

//==============================================================================
bool PianoAudioProcessor::hasEditor() const
{
    return true;
}

juce::AudioProcessorEditor* PianoAudioProcessor::createEditor()
{
    return new PianoAudioProcessorEditor (*this);
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new PianoAudioProcessor();
}
