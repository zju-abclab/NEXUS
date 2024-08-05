#pragma once

#include <NTL/RR.h>
// #include "Remez.h"
#include "MinicompRemez.h"
// #include "func.h"
#include "MinicompFunc.h"
#include "PolyUpdate.h"

using namespace minicomp;

void testRemez();
void testReLU(size_t deg);
RR GetError(int d, RR t, bool is_first_function, int type, RR scale);
// void GetCoeff(int d, RR t, RR *coeff, bool is_first_function, int type, RR scale);
RR GetErrorCoeff(int d, RR t, vector<RR> &coeff, bool is_first_function, int type, RR scale);
// RR GetErrorCoeff_extpoint(int d, RR t, vector<RR> &coeff, vector<RR> &ext_xpos, bool is_first_function, int type, RR scale);
