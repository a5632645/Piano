#include "PluginProcessor.h"
#include "PluginEditor.h"

float defaultPreset[9] = { -8, -8, -6, 0, 0, 0, 0, 0, 0 };

static juce::String onOffTextFunction (const gin::Parameter&, float v)     { return v > 0.5f ? "On" : "Off"; }
static juce::String pVolTextFunction (const gin::Parameter&, float v)      { return v > 0.5f ? "Soft" : "Norm"; }
static juce::String pDecayTextFunction (const gin::Parameter&, float v)    { return v > 0.5f ? "Fast" : "Slow"; }
static juce::String pHarmTextFunction (const gin::Parameter&, float v)     { return v > 0.5f ? "2nd" : "3rd"; }
static juce::String percentTextFunction (const gin::Parameter&, float v)   { return juce::String (juce::roundToInt(v * 100)) + "%"; }

static juce::String vcTextFunction (const gin::Parameter&, float v)
{
    switch (juce::roundToInt (v))
    {
        case 0: return "V1";
        case 1: return "C1";
        case 2: return "V2";
        case 3: return "C2";
        case 4: return "V3";
        case 5: return "C3";
        default: return "";
    }
}

static juce::String lesTextFunction (const gin::Parameter&, float v)
{
    switch (juce::roundToInt (v))
    {
        case 0: return "Stop";
        case 1: return "Slow";
        case 2: return "Fast";
        default: return "";
    }
}

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
