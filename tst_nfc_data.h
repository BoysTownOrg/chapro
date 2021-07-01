// cha_data.h - data size = 4784 bytes
#ifndef CHA_DATA_H
#define CHA_DATA_H

static CHA_DATA cha_magic[4] = {0x55530, 0x68131, 48, 4784};
static CHA_DATA p00[64] = { // _size
    256, 128, 256, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 48, 512, 512, 1024, 1024, 1024};
static CHA_DATA p01[32] = { // _ivar
    32, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 128, 12, 0, 2};
static double p02[32] = { // _dvar
    24};
// empty array: p03
// empty array: p04
// empty array: p05
// empty array: p06
// empty array: p07
// empty array: p08
// empty array: p09
// empty array: p10
// empty array: p11
// empty array: p12
// empty array: p13
// empty array: p14
// empty array: p15
// empty array: p16
// empty array: p17
// empty array: p18
// empty array: p19
// empty array: p20
// empty array: p21
// empty array: p22
// empty array: p23
// empty array: p24
// empty array: p25
// empty array: p26
// empty array: p27
// empty array: p28
// empty array: p29
// empty array: p30
// empty array: p31
// empty array: p32
// empty array: p33
// empty array: p34
// empty array: p35
// empty array: p36
// empty array: p37
// empty array: p38
// empty array: p39
// empty array: p40
// empty array: p41
static CHA_DATA p42[12] = {
    0x00000020, 0x00000025, 0x0000002B, 0x00000031, 0x00000038, 0x0000003F,
    0x00000048, 0x00000051, 0x0000005B, 0x00000066, 0x00000073, 0x00000080};
static CHA_DATA p43[128] = {
    0x3D97B426, 0x3D98C122, 0x3D9BE770, 0x3DA1251E, 0x3DA876F2, 0x3DB1D867,
    0x3DBD43B6, 0x3DCAB1D1, 0x3DDA1A72, 0x3DEB7418, 0x3DFEB40F, 0x3E09E73C,
    0x3E155B27, 0x3E21AEBA, 0x3E2EDA59, 0x3E3CD5E6, 0x3E4B98C0, 0x3E5B19CF,
    0x3E6B4F80, 0x3E7C2FD5, 0x3E86D833, 0x3E8FE335, 0x3E99335A, 0x3EA2C2E6,
    0x3EAC8BF2, 0x3EB68876, 0x3EC0B24A, 0x3ECB0328, 0x3ED574B3, 0x3EE0007C,
    0x3EEAA000, 0x3EF54CB4, 0x3F000000, 0x3F0559A6, 0x3F0AB000, 0x3F0FFFC2,
    0x3F1545A6, 0x3F1A7E6C, 0x3F1FA6DB, 0x3F24BBC5, 0x3F29BA06, 0x3F2E9E8D,
    0x3F336653, 0x3F380E66, 0x3F3C93E7, 0x3F40F40B, 0x3F452C20, 0x3F49398C,
    0x3F4D19D0, 0x3F50CA86, 0x3F54496A, 0x3F579451, 0x3F5AA936, 0x3F5D8631,
    0x3F60297E, 0x3F62917D, 0x3F64BCB2, 0x3F66A9C6, 0x3F685789, 0x3F69C4F3,
    0x3F6AF122, 0x3F6BDB5D, 0x3F6C8312, 0x3F6CE7DC, 0x3F6D097B, 0x3F6CE7DC,
    0x3F6C8312, 0x3F6BDB5D, 0x3F6AF122, 0x3F69C4F3, 0x3F685789, 0x3F66A9C6,
    0x3F64BCB2, 0x3F62917D, 0x3F60297E, 0x3F5D8631, 0x3F5AA936, 0x3F579451,
    0x3F54496A, 0x3F50CA86, 0x3F4D19D0, 0x3F49398C, 0x3F452C20, 0x3F40F40B,
    0x3F3C93E7, 0x3F380E66, 0x3F336653, 0x3F2E9E8D, 0x3F29BA06, 0x3F24BBC5,
    0x3F1FA6DB, 0x3F1A7E6C, 0x3F1545A6, 0x3F0FFFC2, 0x3F0AB000, 0x3F0559A6,
    0x3F000000, 0x3EF54CB4, 0x3EEAA000, 0x3EE0007C, 0x3ED574B3, 0x3ECB0328,
    0x3EC0B24A, 0x3EB68876, 0x3EAC8BF2, 0x3EA2C2E6, 0x3E99335A, 0x3E8FE335,
    0x3E86D833, 0x3E7C2FD5, 0x3E6B4F80, 0x3E5B19CF, 0x3E4B98C0, 0x3E3CD5E6,
    0x3E2EDA59, 0x3E21AEBA, 0x3E155B27, 0x3E09E73C, 0x3DFEB40F, 0x3DEB7418,
    0x3DDA1A72, 0x3DCAB1D1, 0x3DBD43B6, 0x3DB1D867, 0x3DA876F2, 0x3DA1251E,
    0x3D9BE770, 0x3D98C122};
static CHA_DATA p44[128] = {0};
static CHA_DATA p45[256] = {0};
static CHA_DATA p46[256] = {0};
static CHA_DATA p47[256] = {0};

static CHA_DATA *cha_data[64] = {
    (CHA_DATA *)p00, (CHA_DATA *)p01, (CHA_DATA *)p02, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, p42, p43, p44, p45, p46, p47};

#endif // CHA_DATA_H
