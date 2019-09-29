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
#include <cmath>


GrateTime::GrateTime(int year, int month, int day, int hour, int minute, int second) {
    timeinfo.tm_year = year - 1900;
    timeinfo.tm_mon = month;
    timeinfo.tm_mday = day;
    timeinfo.tm_hour = hour;
    timeinfo.tm_min = minute;
    timeinfo.tm_sec = second;
    timeinfo.tm_isdst = 0;
    std::mktime(&timeinfo);
}

// sets the date
void GrateTime::setDate(int year, int month, int day) {
    timeinfo.tm_year = year - 1900;
    timeinfo.tm_mon = month;
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

void GrateTime::setExcelTime(double nSerialDate)
{
    int i,j,l;
    double k,m,n;

       // https://www.codeproject.com/Articles/2750/Excel-Serial-Date-to-Day-Month-Year-and-Vice-Versa
   if (nSerialDate == 60)
   {
          timeinfo.tm_mday = 29;
          timeinfo.tm_mon = 1;
          timeinfo.tm_year = 0;
          std::mktime(&timeinfo);
          return;
   }
   else if (nSerialDate < 60)
   {
          // Because of the 29-02-1900 bug, any serial date
          // under 60 is one off... Compensate.
          nSerialDate++;
   }
   // Modified Julian to DMY calculation with an addition of 2415019

   l = floor(nSerialDate + 68569 + 2415019);
   n = int((4 * l) / 146097);
   l = l - int((146097 * n + 3) / 4);
   i = int((4000 * (l + 1)) / 1461001);
   l = l - int((1461 * i) / 4) + 31;
   j = int((80 * l) / 2447);
   timeinfo.tm_mday = l - int((2447 * j) / 80);
   l = int(j / 11);
   timeinfo.tm_mon = j + 2 - (12 * l) - 1;     // struct tm numbers the months 0 to 11, not 1 to 12
   timeinfo.tm_year = ( 100 * (n - 49) + i + l) - 1900;

   k = (nSerialDate - floor(nSerialDate)) * 24 + 1;
   timeinfo.tm_hour = floor(k);
   m = (k - timeinfo.tm_hour) * 60;
   timeinfo.tm_min = floor(m);
   n = (m - timeinfo.tm_min) * 60;
   timeinfo.tm_sec = round(n);

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
