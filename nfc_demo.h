/*
   NFC_DEMO

   Created: Ira Ray Jenkins, Creare, June 2021

   Based on MyAudioAlgorithm by Chip Audette, OpenAudio April 2017
   Purpose: Demonstrate a port of tst_nfc.c from ChaPro to Tympan.

   MIT License.  use at your own risk.
*/

#include <arm_math.h> //ARM DSP extensions.  https://www.keil.com/pack/doc/CMSIS/DSP/html/index.html
#include "AudioStream_F32.h"

#include <Arduino.h>

#include "chapro.h"

#include "tst_nfc_data.h"

#define MAX_MSG 256

class NFC_DEMO : public AudioStream_F32
{
public:
    //constructor
    NFC_DEMO(const AudioSettings_F32 &settings) : AudioStream_F32(1, inputQueueArray_f32){
                                                      //do any setup activities here

                                                      //if you need the sample rate, it is: fs_Hz = settings.sample_rate_Hz;
                                                      //if you need the block size, it is: n = settings.audio_block_samples;
                                                  };

    //here's the method that is called automatically by the Teensy Audio Library
    void update(void)
    {
        //Serial.println("AudioEffectMine_F32: doing update()");  //for debugging.
        audio_block_f32_t *audio_block;
        audio_block = AudioStream_F32::receiveWritable_f32();
        if (!audio_block)
            return;

        //do your work
        applyMyAlgorithm(audio_block);

        ///transmit the block and release memory
        AudioStream_F32::transmit(audio_block);
        AudioStream_F32::release(audio_block);
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    // Here is where you can add your algorithm.
    // This function gets called block-wise...which is usually hard-coded to every 128 samples
    void applyMyAlgorithm(audio_block_f32_t *audio_block)
    {
        float *x = audio_block->data;
        int cs = audio_block->length;

        cha_nfc_process(this->cp, x, x, cs);

    } //end of applyMyAlgorithms
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    void setup(int chunk_size)
    {
        Serial.println("Setup now...\n");
    }

private:
    //state-related variables
    audio_block_f32_t *inputQueueArray_f32[1]; //memory pointer for the input to this module

    CHA_PTR cp = (CHA_PTR)cha_data;

}; //end class definition for NFC_DEMO
