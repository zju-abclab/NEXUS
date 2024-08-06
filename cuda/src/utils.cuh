#pragma once

#include <chrono>
#include <iostream>

namespace nexus {
using namespace std;
using namespace std::chrono;

class Timer {
 private:
  std::chrono::_V2::system_clock::time_point start_;
  std::chrono::_V2::system_clock::time_point end_;

  template <typename T = milliseconds>
  inline void print_duration(const char *text) {
    cout << text << duration_cast<T>(end_ - start_).count() << endl;
  }

 public:
  Timer() {
    start_ = high_resolution_clock::now();
  }

  inline void start() {
    start_ = high_resolution_clock::now();
  }

  template <typename T = milliseconds>
  inline void stop(const char *text) {
    end_ = high_resolution_clock::now();
    print_duration<T>(text);
  }

  inline void stop() {
    end_ = high_resolution_clock::now();
  }

  template <typename T = milliseconds>
  inline long duration() {
    return duration_cast<T>(end_ - start_).count() / 1.0;
  }
};
}  // namespace nexus
