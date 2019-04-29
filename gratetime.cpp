/*******************
 *
 *
 *  GRATE 8
 *
 *  Custom Date Time object
 *
 *
 *
*********************/

#include "gratetime.h"
#include <iostream>
#include <ctime>


GrateTime::GrateTime(int year, int month, int day, int hour, int minute, int second) {
    setDate(year, month, day);
    setTime(hour, minute, second);
}

// sets the date
void GrateTime::setDate(int year, int month, int day) {
    timeinfo.tm_year = year - 1900;
    timeinfo.tm_mon = month - 1;
    timeinfo.tm_mday = day;

    // normalise it
    std::mktime(&timeinfo);
}

// sets the time
void GrateTime::setTime(int hour, int minute, int second) {
    timeinfo.tm_hour = hour;
    timeinfo.tm_min = minute;
    timeinfo.tm_sec = second;

    // normalise it
    std::mktime(&timeinfo);
}

// prints a string representation to stdout
void GrateTime::print() {
    std::cout << std::asctime(&timeinfo) << std::endl;
}

// adds the given number of seconds in-place
void GrateTime::addSecs(int secondsToAdd) {
    timeinfo.tm_sec += secondsToAdd;

    // normalise it
    std::mktime(&timeinfo);
}

// returns the date-time as time_t type
std::time_t GrateTime::getTime_t() {
    return mktime(&timeinfo);
}

// calculates the number of seconds to the given GrateTime object
int GrateTime::secsTo(GrateTime futureTime) {
    return std::difftime(futureTime.getTime_t(), getTime_t());
}
