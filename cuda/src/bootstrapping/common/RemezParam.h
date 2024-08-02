#pragma once

#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

class RemezParam {
 public:
  double log_scan_step_diff = 9.5;
  long binary_prec = 10;
  long RR_prec = 1000;
  long log_approx_degree = 120;
  long log_round_prec = 100;

  RemezParam() {}
  RemezParam(double scd, long bi, long rr, long app, long round)
      : log_scan_step_diff(scd), binary_prec(bi), RR_prec(rr), log_approx_degree(app), log_round_prec(round) {}

  RemezParam(const RemezParam& copy)
      : log_scan_step_diff(copy.log_scan_step_diff), binary_prec(copy.binary_prec), RR_prec(copy.RR_prec), log_approx_degree(copy.log_approx_degree), log_round_prec(copy.log_round_prec) {}
};
