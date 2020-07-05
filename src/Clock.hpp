#ifndef _CLOCK_HPP
#define _CLOCK_HPP

#include <iostream>
#include <iomanip>
#include <chrono>
#include <math.h>
class Clock {
protected:
  std::chrono::steady_clock::time_point startTime;
public:
  Clock() {
    tick();
  }
  void tick() {
    startTime = std::chrono::steady_clock::now();
  }
  float tock() {
    std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();;
    return float(std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()) / 1e6f;
  }
  // Print elapsed time with newline
  void ptock() {
    float elapsed  = tock();
    double hours   = floor(elapsed/(60*60));
    double minutes = floor((elapsed - hours*60*60)/60.);
    double seconds = elapsed - hours*60*60 - minutes*60;
    if (hours > 1)
      std::cout << "TOOK " << hours << " HOURS " << minutes << " MINUTES, " << seconds << " SECONDS" << std::endl;
    else if (minutes > 1)
      std::cout << "TOOK " << minutes << " MINUTES, " << seconds << " SECONDS" << std::endl;
    else
      std::cout << "TOOK " << elapsed << " SECONDS" << std::endl;
  }
};

#endif
