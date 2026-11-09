
#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <string>
#include <ctime>

/* Timer is for measuring wall clock time and CPU time of different code sections. Stop is used to make a checkpoint in time and reset for restarting the timer. */
class Timer
{

public:
    Timer();
    void reset();
    void stop();

    std::string getCPUTimeStr();
    std::string getWallTimeStr();
    std::string getStartTimeStr();

private:
    time_t wallEndTime, wallStartTime;
    clock_t cpuEndTime, cpuStartTime;
    std::string getTimeStr(double time, int numDecimals);
};

/* StageTimer reports the wall clock and CPU time spent in its scope to stderr
   when the verbosity is above 2. Used to see which processing stage dominates
   on large inputs. */
class StageTimer
{

public:
    explicit StageTimer(const char* name);
    ~StageTimer();

private:
    const char* name_;
    std::chrono::steady_clock::time_point wallStart_;
    clock_t cpuStart_;
};

#endif

