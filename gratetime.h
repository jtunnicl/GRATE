#ifndef _GRATETIME_H_
#define _GRATETIME_H_

#include <ctime>


// Date Time object so we don't have to rely on QDateTime
class GrateTime
{
    public:
        GrateTime(int year = 2000, int month = 1, int day = 1, int hour = 0, int minute = 0, int second = 0);

        void addSecs(int secondsToAdd);  // add seconds onto current time
        int secsTo(GrateTime futureTime);  // calculate number of seconds to the given GrateTime object
        std::time_t getTime_t();  // returns the time in time_t format
        void setDate(int year, int month, int day);  // set the date
        void setTime(int hour, int minute, int second);  // set the time
        void setExcelTime(double nSerialDate);  // Convert from Excel serial date
        void print();  // print the current date time

    private:
        std::tm timeinfo;
};

#endif
