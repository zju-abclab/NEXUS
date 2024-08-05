#include "RemezApp.h"

// d is odd!!
void testReLU(size_t deg) {
  RR maxerr;
  //	long prec = 300;
  long prec = 50;
  long iteration = 200;
  long RRprec;
  int num;
  RR sc;
  //	long d = 300;
  long type = 1;
  RR scale = RR(1);

  size_t inter_num = 1;
  vector<RR> inter_start, inter_end;
  vector<RR> coeff, chebcoeff;
  inter_start.emplace_back(RR(-1));
  inter_end.emplace_back(RR(1));

  num = 1000;
  //	num = 2000;
  RRprec = 300;
  //	RRprec = 2000;
  RR::SetPrecision(RRprec);
  sc = RR(2.0) / RR(num);

  Remez rm(ReLU, inter_num, inter_start, inter_end, sc, prec, deg, iteration, type, scale, true);
  rm.test();
  maxerr = rm.getmaxerr();
  cout << "maxerr: " << maxerr << endl;

  rm.getcoeff(chebcoeff);
  Polynomial poly(deg, chebcoeff, "cheb");
  poly.cheb_to_power();
  poly.get_coeff(coeff);
  //	cout << "coefficients: " << endl;
  ofstream out1("../result/ReLUpower.txt"), out2("../result/ReLUcheb.txt");
  for (auto num : coeff) out1 << num << endl;
  for (auto num : chebcoeff) out2 << num << endl;
  cout << endl;

  // test for coefficients
  RR max = RR(0);
  for (RR x = RR(-1); x < RR(1); x += RR(0.0001))
    if (abs(poly.evaluate(x) - ReLU(x)) > max) max = abs(poly.evaluate(x) - ReLU(x));
  cout << "power basis approx maxerr: " << max << endl;
  max = RR(0);
  for (RR x = RR(-1); x < RR(1); x += RR(0.0001))
    if (abs(poly.evaluate_cheb(x) - ReLU(x)) > max) max = abs(poly.evaluate_cheb(x) - ReLU(x));
  cout << "cheb basis approx maxerr: " << max << endl;
}
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
/*
void GetCoeff(int d, RR t, RR coeff[], bool is_first_function, int type, RR scale)
{
        Remez *rm;
        RR maxerr;
        long prec = 40;
        long iteration = 20;
        long RRprec;
        int num;
        RR sc;

        //  num & RRprec setting!!!
        num = 2000;
        RRprec = 800;

        RR::SetPrecision(RRprec);
        sc = t/static_cast<RR>(num);
        rm = new Remez(t, is_first_function, sc, prec, d+1, iteration, type, scale);
        rm->test();
        rm->getcoeff(coeff);

        delete rm;
}
*/
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
/*
RR GetErrorCoeff_extpoint(int d, RR t, vector<RR> &coeff, vector<RR> &ext_xpos, bool is_first_function, int type, RR scale)
{
        Remez *rm;
        RR maxerr;
        long prec = 40;
        long iteration = 20;
        long RRprec;
        int num;;
        RR sc;

        //  num & RRprec setting!!!
        num = 1000;
        RRprec = 300;

        RR::SetPrecision(RRprec);
        sc = t/static_cast<RR>(num);
        rm = new Remez(t, is_first_function, sc, prec, d+1, iteration, type, scale);
        rm->test();
        rm->getcoeff(coeff);
        rm->getext_xpos(ext_xpos);
        maxerr = rm->getmaxerr();

        delete rm;
        return maxerr;
}
*/
