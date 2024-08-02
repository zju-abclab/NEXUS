#pragma once

#include "common/Remez.cuh"

class RemezCos : public boot::Remez {
 public:
  long scale_factor;

  RemezCos(RemezParam _param, long _boundary_K, long _log_width, long _deg, long _scale_factor)
      : boot::Remez(_param, _boundary_K, _log_width, _deg), scale_factor(_scale_factor) {}

  RR function_value(RR x) {
    if (scale_factor % 2 == 0)
      return cos(2 * ComputePi_RR() * (x - 0.25) / scale_factor);
    else
      return sin(2 * ComputePi_RR() * x / scale_factor);
  }
};
