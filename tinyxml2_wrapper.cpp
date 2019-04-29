
#include "tinyxml2_wrapper.h"
#include <sstream>
#include <iostream>

using namespace tinyxml2;

double getDoubleValue(XMLElement *e, const char *name) {
    // get the child element
    XMLElement *child = e->FirstChildElement(name);
    if (child == NULL) {
        std::stringstream error_stream;
        error_stream << "Error getting child element: " << name;
        std::string msg = error_stream.str();
        throw msg;
    }

    // get the value
    double value;
    if (child->QueryDoubleText(&value) != XML_SUCCESS) {
        std::stringstream error_stream;
        error_stream << "Error getting double value for child element: " << name;
        std::string msg = error_stream.str();
        throw msg;
    }

    return value;
}

float getFloatValue(XMLElement *e, const char *name) {
    // get the child element
    XMLElement *child = e->FirstChildElement(name);
    if (child == NULL) {
        std::stringstream error_stream;
        error_stream << "Error getting child element: " << name;
        std::string msg = error_stream.str();
        throw msg;
    }

    // get the value
    float value;
    if (child->QueryFloatText(&value) != XML_SUCCESS) {
        std::stringstream error_stream;
        error_stream << "Error getting double value for child element: " << name;
        std::string msg = error_stream.str();
        throw msg;
    }

    return value;
}

int getIntValue(XMLElement *e, const char *name) {
    // get the child element
    XMLElement *child = e->FirstChildElement(name);
    if (child == NULL) {
        std::stringstream error_stream;
        error_stream << "Error getting child element: " << name;
        std::string msg = error_stream.str();
        throw msg;
    }

    // get the value
    int value;
    if (child->QueryIntText(&value) != XML_SUCCESS) {
        std::stringstream error_stream;
        error_stream << "Error getting int value for child element: " << name;
        std::string msg = error_stream.str();
        throw msg;
    }

    return value;
}
