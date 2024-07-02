#pragma once

#include <NTL/mat_RR.h>

#include <string>

#include "Choosemax.cuh"
#include "MinicompFunc.cuh"
#include "Point.cuh"
#define _USE_MATH_DEFINES

using namespace std;
using namespace NTL;

namespace minicomp {
class Remez {
 private:
  // input variables
  long deg, iter, prec;
  int type;
  RR scale, sc;
  //	RR sc, *coeff, *temp_ext, maxerr;
  vector<RR> coeff;
  RR(*func)
  (RR);
  size_t inter_num;
  vector<RR> inter_start, inter_end;
  //	bool is_first_function;
  bool is_opt_sampling;

  // variables for algorithm
  Point *sam, *ext;
  RR maxerr;
  long ext_count;
  int *ext_maxsum_index;
  RR *temp_ext;
  //	string filename;
  vec_RR v, w, v_0;
  mat_RR m;

 public:
  Remez(RR (*_func)(RR), size_t _inter_num, vector<RR> _inter_start, vector<RR> _inter_end, RR _sc, long _prec, long _deg, long _iter, int _type, RR _scale, bool _is_opt_sampling);
  ~Remez();
  void initialize();
  void getcoeffwitherr();
  int getextreme();
  void choosemaxs();
  void printgraph(ofstream &out, RR start, RR end, RR sc);
  bool test();
  RR getmaxerr();
  void getcoeff(RR _coeff[]);
  void getcoeff(vector<RR> &_coeff);
  void getext_xpos(vector<RR> &_ext_xpos);
};
}  // namespace minicomp
