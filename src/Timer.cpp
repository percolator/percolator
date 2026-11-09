

#include "Timer.h"

#include <iostream>
#include <sstream>

#include "Globals.h"

Timer::Timer(){
    reset();
}

void Timer::reset()
{
    cpuStartTime = cpuEndTime = clock();
    time(&wallStartTime);
    time(&wallEndTime);
}
void Timer::stop()
{
    cpuEndTime = clock();
    time(&wallEndTime);
}

std::string Timer::getCPUTimeStr()
{
    double cpuTime = ((double)(cpuEndTime - cpuStartTime)) / (double)CLOCKS_PER_SEC;
    return getTimeStr(cpuTime, 4);
}

std::string Timer::getWallTimeStr()
{
    return getTimeStr(difftime(wallEndTime, wallStartTime), 0);
}

std::string Timer::getTimeStr(double time, int numDecimals)
{
    std::stringstream timeStr;
    timeStr.precision(numDecimals);
    timeStr << std::fixed << time;
    return timeStr.str();
}

std::string Timer::getStartTimeStr(){
    return ctime(&wallStartTime);
}

StageTimer::StageTimer(const char* name)
    : name_(name), wallStart_(std::chrono::steady_clock::now()), cpuStart_(clock()) {}

StageTimer::~StageTimer()
{
    if (VERB > 2) {
        double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wallStart_).count();
        double cpu = ((double)(clock() - cpuStart_)) / (double)CLOCKS_PER_SEC;
        std::cerr << "Stage " << name_ << ": " << wall << " s wall, " << cpu << " s cpu" << std::endl;
    }
}
