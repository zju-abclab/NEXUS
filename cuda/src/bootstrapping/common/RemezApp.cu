#include "RemezApp.cuh"

RR GetError(int d, RR t, bool is_first_function, int type, RR scale) {
  Remez *rm;
  RR maxerr;
  long prec = 40;
  long iteration = 20;
  long RRprec;
  int num;
  RR sc;

  size_t inter_num = 2;
  vector<RR> inter_start, inter_end;
  if (is_first_function == true) {
    inter_start.emplace_back(RR(-1));
    inter_start.emplace_back(RR(1) - t);
    inter_end.emplace_back(RR(-1) + t);
    inter_end.emplace_back(RR(1));
  } else {
    inter_start.emplace_back(RR(-1) - t);
    inter_start.emplace_back(RR(1) - t);
    inter_end.emplace_back(RR(-1) + t);
    inter_end.emplace_back(RR(1) + t);
  }

  //  num & RRprec setting!!!
  num = 500;
  RRprec = 1500;  // for deg 64. RRprec >= 1500 approximately. But, if tau > 1/2, RRprec = 300 is ok.

  RR::SetPrecision(RRprec);
  sc = t / static_cast<RR>(num);
  //	rm = new Remez(t, is_first_function, sc, prec, d+1, iteration, type, scale);
  rm = new Remez(sgn, inter_num, inter_start, inter_end, sc, prec, d + 1, iteration, type, scale, true);
  rm->test();
  maxerr = rm->getmaxerr();

  delete rm;
  return maxerr;
}

RR GetErrorCoeff(int d, RR t, vector<RR> &coeff, bool is_first_function, int type, RR scale) {
  Remez *rm;
  RR maxerr;
  long prec = 40;
  long iteration = 25;
  long RRprec;
  int num;
  ;
  RR sc;

  size_t inter_num = 2;
  vector<RR> inter_start, inter_end;
  if (is_first_function == true) {
    inter_start.emplace_back(RR(-1));
    inter_start.emplace_back(RR(1) - t);
    inter_end.emplace_back(RR(-1) + t);
    inter_end.emplace_back(RR(1));
  } else {
    inter_start.emplace_back(RR(-1) - t);
    inter_start.emplace_back(RR(1) - t);
    inter_end.emplace_back(RR(-1) + t);
    inter_end.emplace_back(RR(1) + t);
  }

  //  num & RRprec setting!!!
  num = 1000;
  RRprec = 300;

  RR::SetPrecision(RRprec);
  sc = t / static_cast<RR>(num);
  //	rm = new Remez(t, is_first_function, sc, prec, d+1, iteration, type, scale);
  rm = new Remez(sgn, inter_num, inter_start, inter_end, sc, prec, d + 1, iteration, type, scale, true);
  rm->test();
  rm->getcoeff(coeff);
  maxerr = rm->getmaxerr();

  delete rm;
  return maxerr;
}
