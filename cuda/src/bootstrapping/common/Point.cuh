#pragma once
#include <NTL/RR.h>

#include <iostream>

using namespace std;
using namespace NTL;

class Point {
 public:
  RR x;
  RR y;

  long locmm;

  Point() = default;
  Point(RR _x, RR _y);
};
