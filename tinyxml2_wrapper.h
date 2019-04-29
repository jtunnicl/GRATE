#ifndef _TINYXML2_WRAPPER_H_
#define _TINYXML2_WRAPPER_H_

#include "tinyxml2/tinyxml2.h"

using namespace tinyxml2;

double getDoubleValue(XMLElement *e, const char *name);
float getFloatValue(XMLElement *e, const char *name);
int getIntValue(XMLElement *e, const char *name);

#endif
