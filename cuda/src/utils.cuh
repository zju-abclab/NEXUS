#pragma once
#include <chrono>
#include <iostream>

using namespace std;
using namespace std::chrono;

namespace nexus {
class Timer {
 private:
  std::chrono::_V2::system_clock::time_point start;
  std::chrono::_V2::system_clock::time_point end;

  template <typename T = milliseconds>
  inline void print_duration(const char *text) {
    cout << text << duration_cast<T>(end - start).count() << endl;
  }

 public:
  Timer() {
    start = high_resolution_clock::now();
  }

  template <typename T = milliseconds>
  inline void stop(const char *text) {
    end = high_resolution_clock::now();
    print_duration<T>(text);
  }

  inline void stop() {
    end = high_resolution_clock::now();
  }

  template <typename T = milliseconds>
  inline long duration() {
    return duration_cast<T>(end - start).count() / 1.0;
  }
};
}  // namespace nexus
