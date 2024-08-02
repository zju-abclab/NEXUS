#pragma once

#include "common/Remez.cuh"

class RemezArcsin : public boot::Remez {
 public:
  RemezArcsin(RemezParam _param, double _log_width, long _deg)
      : boot::Remez(_param, 1, _log_width, _deg) {}

  RR function_value(RR x) {
    return 1 / (2 * ComputePi_RR()) * arcsin(x);
  }
};
