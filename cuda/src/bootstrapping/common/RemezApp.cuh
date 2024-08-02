#pragma once

#include <NTL/RR.h>

#include "MinicompFunc.cuh"
#include "MinicompRemez.cuh"
#include "PolyUpdate.cuh"

using namespace minicomp;

RR GetError(int d, RR t, bool is_first_function, int type, RR scale);
RR GetErrorCoeff(int d, RR t, vector<RR> &coeff, bool is_first_function, int type, RR scale);
