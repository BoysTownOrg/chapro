#include <Arduino.h>
extern "C"
{
#include "seriallog.h"
}

void seriallog(const char *msg)
{
    Serial.println(msg);
}
