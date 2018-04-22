// cha_data.h - array size = 5840 bytes
#ifndef CHA_DATA_H
#define CHA_DATA_H

static CHA_DATA p00[      64] = { // _size
               256,        64,       128,      2048,       160,       160,       160,      1888,
                32,        32,         0,         0,         0,         0,         0,         0,
                 0,       512,       400
};
static CHA_DATA p01[      16] = { // _ivar
                32,         8,         5,        59,       128,         0,         0,       100
};
static double   p02[      16] = { // _dvar
                      0,              0,          24000,              0,              0,
                      0,              0,              0,              0,              0,
                  0.001,          0.984
};
static CHA_DATA p03[     512] = {         0};
static CHA_DATA p04[      40] = {
        0x3620D502,0x3720D502,0x37713F84,0x3720D502,0x3620D502,0x39F5BEF2,0x00000000,0xBA75BEF2,
        0x00000000,0x39F5BEF2,0x3AA38234,0x00000000,0xBB238234,0x00000000,0x3AA38234,0x3B513956,
        0x00000000,0xBBD13956,0x00000000,0x3B513956,0x3C0124F2,0x00000000,0xBC8124F2,0x00000000,
        0x3C0124F2,0x3C99FF26,0x00000000,0xBD19FF26,0x00000000,0x3C99FF26,0x3D308AD6,0x00000000,
        0xBDB08AD6,0x00000000,0x3D308AD6,0x3E18B5CF,0xBF18B5CF,0x3F6510B6,0xBF18B5CF,0x3E18B5CF
};
static CHA_DATA p05[      40] = {
        0x3F800000,0xC0727F33,0x40AC72CE,0xC05A380A,0x3F4F490A,0x3F800000,0xC07AA159,0x40B8B9DE,
        0xC072F4A5,0x3F7090C7,0x3F800000,0xC07600D5,0x40B2FF21,0xC069B631,0x3F671295,0x3F800000,
        0xC06CE446,0x40A87740,0xC059F28B,0x3F58C157,0x3F800000,0xC05A5B66,0x4094F801,0xC03EE2F0,
        0x3F43E08B,0x3F800000,0xC0343075,0x4064EC82,0xC01115DC,0x3F26D0EE,0x3F800000,0xBFD056F9,
        0x4000A17F,0xBF928031,0x3F015F3D,0x3F800000,0xBF1DC13F,0x3F1B35F6,0xBE0DA672,0x3CCECD2A
};
static CHA_DATA p06[      40] = {         0};
static CHA_DATA p07[     472] = {         0};
static CHA_DATA p08[       8] = {
        0x3F800000,0x3F034AFA,0x3F3BC5A1,0x3F453F33,0x3F2DE439,0x3F47FF0A,0x3F5F7AF5,0xBF800000
};
static CHA_DATA p09[       8] = {
        0x41B80000,0x3F800000,0x41B00000,0x42100000,0x42300000,0x42480000,0x42540000,0x42680000
};
// empty array ->     p10
// empty array ->     p11
// empty array ->     p12
// empty array ->     p13
// empty array ->     p14
// empty array ->     p15
// empty array ->     p16
static CHA_DATA p17[     128] = {         0};
static CHA_DATA p18[     100] = {         0};

static CHA_DATA *cha_data[NPTR] = {
    (CHA_DATA *)p00,(CHA_DATA *)p01,(CHA_DATA *)p02,
     p03, p04, p05, p06, p07, p08, p09,NULL,NULL,NULL,NULL,NULL,NULL,NULL, p17, p18
};

#endif // CHA_DATA_H
