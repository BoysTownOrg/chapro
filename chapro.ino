#include <Arduino.h>
#include <Tympan_Library.h>
#include "AudioSDPlayer_F32.h"
#include "NFC_DEMO.h"
#include <SD.h>

#include "chapro.h"

//set the sample rate and block size
const float sample_rate_Hz = (int)(44100); //this value puts 37000Hz right at a FFT bin (or other frequencies in the table in AudioOutputI2S_F32)
const int audio_block_samples = 128;       //do not make bigger than AUDIO_BLOCK_SAMPLES from AudioStream.h (which is 128)  Must be 128 for SD recording.
AudioSettings_F32 audio_settings(sample_rate_Hz, audio_block_samples);

Tympan myTympan(TympanRev::E, audio_settings);

//create audio objects
AudioSDPlayer_F32 audioSDPlayer(audio_settings);
NFC_DEMO effect1(audio_settings);
NFC_DEMO effect2(audio_settings);
AudioOutputI2S_F32 audioOutput(audio_settings);

//create audio connections
AudioConnection_F32 patchCord1(audioSDPlayer, 0, effect1, 0);
AudioConnection_F32 patchCord2(audioSDPlayer, 1, effect2, 0);
AudioConnection_F32 patchCord3(effect1, 0, audioOutput, 0);
AudioConnection_F32 patchCord4(effect2, 0, audioOutput, 1);

unsigned long end_millis = 0;
String filename = "CAT.WAV"; // filenames are always uppercase 8.3 format

void setup()
{
    myTympan.beginBothSerial();
    delay(5000);
    myTympan.print("SDWavPlayer");
    myTympan.println(": setup():...");
    myTympan.print("Sample Rate (Hz): ");
    myTympan.println(audio_settings.sample_rate_Hz);
    myTympan.print("Audio Block Size (samples): ");
    myTympan.println(audio_settings.audio_block_samples);

    // Audio connections require memory to work.
    AudioMemory_F32(20, audio_settings);

    effect1.setup(audio_block_samples);

    // Start the Tympan
    myTympan.enable();
    myTympan.volume(0.5);

    //prepare SD player
    audioSDPlayer.begin();

    //finish setup
    delay(2000); //stall a second
    Serial.println("Setup complete.");
}

void loop()
{
    //service the audio player
    if (!audioSDPlayer.isPlaying())
    { //wait until previous play is done
        AudioProcessorUsageMaxReset();
        //start playing audio
        myTympan.print("Starting audio player: ");
        myTympan.println(filename);
        audioSDPlayer.play(filename);
    }

    //do other things here, if desired...like, maybe check the volume knob?

    // myTympan.printCPUandMemory(millis(), 1000);

    //stall, just to be nice?
    delay(5000);
}
