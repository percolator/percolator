
#ifndef TIMER_H_
#define TIMER_H_

#include <string>
#include <ctime>

class Timer
{

public:
    Timer();
    void reset();
    void stop();

    std::string getCPUTimeStr();
    std::string getWallTimeStr();
    char* getStartTimeStr();

private:
    time_t wallEndTime, wallStartTime;
    clock_t cpuEndTime, cpuStartTime;
    std::string getTimeStr(double time, int numDecimals);
};

#endif

