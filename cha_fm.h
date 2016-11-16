// cha_fm.h - function prototypes and defines for feedback management
#ifndef CHA_FM_H
#define CHA_FM_H

/*****************************************************/

// feedback module

FUNC(int) cha_feedback_prepare(CHA_PTR, int);
FUNC(void) cha_feedback_manage(CHA_PTR, float *, float *, int);
FUNC(void) cha_feedback_record(CHA_PTR, float *, int);

/*****************************************************/

#define _offset   _reserve

// pointer indices

#define _cc       _offset+0

// integer variable indices

#define _cs       0 

// double variable indices

#define _fs       0

#endif /* CHA_FM_H */
