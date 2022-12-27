
#ifndef PIANO_H
#define PIANO_H

#include <JuceHeader.h>

#define PIANO_MIN_NOTE 21
#define PIANO_MAX_NOTE 108

#undef FDN_REVERB
#define FDN_REVERB 1

const int NUM_NOTES = 128;    // 128 midi notes

#include "filter.h"
#include "dwgs.h"
#include "reverb.h"
#include "hammer.h"

enum {
  MaxDecimation = 16,
  MaxLongDelay = 16,
  TranBufferSize = 32,
  NumParams = 34
};

class Param {
public:
  int tag;
  char name[64];
  char label[64];
};

class Value {
public:
  operator float() const {
    return v;
  }
  float f;
  float v;
};

enum parameters {
  pYoungsModulus = 0,
  pStringDensity,
  pHammerMass,
  pStringTension,
  pStringLength,
  pStringRadius,
  pHammerCompliance,
  pHammerSpringConstant,
  pHammerHysteresis,
  pBridgeImpedance,
  pBridgeHorizontalImpedance,
  pVerticalHorizontalImpedance,
  pHammerPosition,
  pSoundboardSize,
  pStringDecay,
  pStringLopass,
  pDampedStringDecay,
  pDampedStringLopass,
  pSoundboardDecay,
  pSoundboardLopass,
  pLongitudinalGamma,
  pLongitudinalGammaQuadratic,
  pLongitudinalGammaDamped,
  pLongitudinalGammaQuadraticDamped,
  pLongitudinalMix,
  pLongitudinalTransverseMix,
  pVolume,
  pMaxVelocity,
  pStringDetuning,
  pBridgeMass,
  pBridgeSpring,
  pDwgs4,
  pDownsample,
  pLongModes
};

int getParameterIndex(const char *key);

class Piano;

class PianoNote {
public:
	PianoNote(int note, int Fs, Piano *piano);
	~PianoNote();
	float go();
	vec4 go4();
	float goUp();
	float goUpDelayed();
	float goDown();
	float goDownDelayed();
	void triggerOn(float velocity, float *tune);
	void triggerOff();
	bool isActive();
	void deActivate();

	PianoNote *next;
	PianoNote *prev;	
	int note;
  Piano *piano;
  
	bool isDone();
	float maxEnergy;
	float energy;
	Delay<65536> outputdelay;
  static void fillFrequencyTable();
  static double freqTable[NUM_NOTES];    

protected:

  MSD2Filter bridge;


	bool bActive;
	int Fs;
    BiquadHP longHP;
    float alphasb;
	float Z2i;
	float nstringsi;
    float tranBridgeForce;
    float longBridgeForce;
    float longTranBridgeForce;
    float Z;
    float Zhv;
    float vTran;

    int downsample;
    float longForces[MaxDecimation + 1];
    float longForcesK[MaxDecimation + 1];

    float tranForces[TranBufferSize]  __attribute__((aligned(32)));
    float tranForcesH[TranBufferSize]  __attribute__((aligned(32)));
    float longTranForces[TranBufferSize]  __attribute__((aligned(32)));

    float T;
    float mu;
	float f;
    int nstrings;

	dwgs *stringT[3];	
	dwgs *stringHT[3];	
    Hammer *hammer;

    float outUp[8] __attribute__((aligned(32)));
    float outDown[8] __attribute__((aligned(32)));
    int tUp;
    int tDown;
    int upSampleDelayNeeded;
    int downSampleDelayNeeded;

    bool bInit4;
    int tTran;
    int tLong;
    int tTranRead;
    int longDelay;
    int upsample;

    ResampleFIR downSampleFilter;
    ResampleFIR upSampleFilter;
};

class Piano
{
public:
	Piano (int parameters);
	~Piano();

	void addVoice(PianoNote *v);	
	void removeVoice(PianoNote *v);
	void process(float *out, int samples);	
	void process(float **in, float **out, int frameSamples, int offset);
	void process(float **in, float **out, int frameSamples, juce::MidiBuffer&);
    void init(float Fs, int blockSize);
    void triggerOn(int note, float velocity, float *tune);

	//vst crap
	void setParameter (int32_t index, float value);
	float getParameter (int32_t index);

protected:     
    friend class PianoNote;
    Value vals[NumParams];
	PianoNote *voiceList;
	PianoNote *noteArray[NUM_NOTES];
    int blockSize;
	float Fs;
    //int blockSize;
    float *input;

#ifdef FDN_REVERB
  Reverb *soundboard;
#else
  ConvolveReverb<revSize> *soundboard;
#endif
};

void qianoNote(int note, int N, float velocity, double *x, double *tune);
extern "C" void qianoNote(int note, int N, float velocity, float *x, float *tune);
	
#endif
