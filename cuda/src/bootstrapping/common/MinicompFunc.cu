// #include "func.h"
#include "MinicompFunc.cuh"

namespace minicomp {
// [0, q-1) -> (-q/2,q/2]
ZZ mod(ZZ a, ZZ q) {
  if (a > q / 2)
    return a - q;
  else
    return a;
}
long pmod(long a, long b) {
  return (a % b + b) % b;
}
long pow2(long n) {
  long prod = 1;
  for (int i = 0; i < n; i++) prod *= 2;

  return prod;
}
// ceil(x/y). y>0
size_t ceil_divide(size_t x, size_t y) {
  if (y <= 0) std::invalid_argument("y <= 0");
  return (x + (y - 1)) / y;
}
size_t ceil_to_int(double x) {
  return static_cast<size_t>(ceil(x) + 0.5);
}
int floor_to_int(double x) {
  return static_cast<int>(floor(x) + 0.5);
}
// d=-1 implies that n is not power of 2
long log2_long(long n) {
  if (n > 65536 || n <= 0) throw std::out_of_range("n is too large.");
  int d = -1;
  for (int i = 0; i <= 16; i++)
    if (pow2(i) == n) {
      d = i;
      break;
    }

  return d;
}
long num_one(long n) {
  int i = 0, num = 0;
  while (1) {
    if (pow2(i) > n) break;
    i++;
  }
  for (int j = 0; j < i; j++) {
    if (n % 2 == 1) num++;
    n /= 2;
  }
  return num;
}
RR GetApproxError(int d, RR t, bool is_first_function) {
  if (is_first_function == false) {
    return exptoreal(expmaxerr(d, realtoexp(t)));
  } else {
    return exptoreal(expmaxerr(d, realtoexp(t / (2.0 - t))));
  }
}

// only not first function
RR GetInvApproxError(int d, RR t, RR *X, RR *Y, long num) {
  return exptoreal(invexpmaxerr(d, realtoexp(t), X, Y, num));
}
RR GetInvApproxError(int d, RR t, vector<RR> &X, vector<RR> &Y, long num) {
  if (t <= RR(0) || t >= RR(1)) throw std::out_of_range("t should be in (0,1)");
  return exptoreal(invexpmaxerr(d, realtoexp(t), X, Y, num));
}
int enough_mult(int a) {
  //	int m[30] = {0,0,24,24,24,32,32,48,48,48,56,56,64,64,64,72,72,80,80,80,80,80,80,80,80,80,80,80,80,80};
  int m[30];
  for (int i = 0; i < 30; i++) m[i] = 70;
  return m[a];
}
int enough_dep(int a) {
  //	int d[30] = {0,0,12,12,12,16,16,24,24,24,28,28,32,32,32,36,36,40,40,40,44,50,50,50,50,50,50,50,50,50};
  int d[30];
  for (int i = 0; i < 30; i++) d[i] = 40;
  return d[a];
}
int dep(list<int> &lt) {
  list<int>::iterator iter;
  int d[16] = {0, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5};
  int tot_depth = 0;
  for (iter = lt.begin(); iter != lt.end(); iter++) tot_depth += d[((*iter) - 1) / 2];
  return tot_depth;
}
int dep(int deg)  // only odd degrees
{
  int d[32] = {0, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};  // deg: 1,3,5,7,...,63
  return d[(deg - 1) / 2];
}
int mult(int deg)  // only odd degrees
{
  //	int m[16]={0,2,3,4,4,5,6,7,7,8,8,8,10,10,10,10};
  int m[32] = {0, 2, 3, 5, 5, 6, 7, 8, 8, 8, 9, 9, 10, 10, 11, 12, 11, 11, 11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 17};  // deg: 1,3,5,7,...,63
  return m[(deg - 1) / 2];
}

RR exptoreal(RR x) {
  if (x < 0)
    return pow(static_cast<RR>(2.0), x);
  else
    return 1.0 - pow(static_cast<RR>(2.0), -x);
}

RR realtoexp(RR x) {
  if (x <= 0) {
    cout << "error occur" << endl;
    cout << x << endl;
  }
  if (x >= 1)
    return RR(10000);
  else if (x < 0.5)
    return log(x) / log(static_cast<RR>(2));
  else
    return -log(static_cast<RR>(1) - x) / log(static_cast<RR>(2));
}
RR sgn(RR x) {
  if (x > 0)
    return static_cast<RR>(1);
  else if (x < 0)
    return static_cast<RR>(-1);
  else
    return static_cast<RR>(0);
}
RR comp(RR a, RR b) {
  if (a > b)
    return RR(1);
  else if (a < b)
    return RR(0);
  else
    return RR(0.5);
}
RR ReLU(RR x) {
  if (x > RR(0))
    return x;
  else
    return RR(0);
}
RR fracpart(RR x) {
  return x - round(x);
}
double comp(double a, double b) {
  if (a > b)
    return 1.0;
  else if (a < b)
    return 0.0;
  else
    return 0.5;
}
double ReLU(double x) {
  if (x > 0)
    return x;
  else
    return 0.0;
}
// redundant coefficients are not important.
// coeff: power basis coeff or cheb coeff
RR eval(long deg, RR *coeff, RR val, int type, RR scale) {
  RR rtn;
  if (type == 0) {
    RR tmp;
    tmp = 1;
    rtn = coeff[0] * tmp;
    for (int i = 1; i <= deg; i++) {
      tmp = val * tmp;
      rtn += coeff[i] * tmp;
    }
  } else if (type == 1) {
    if (deg == 0) {
      return RR(1);
    } else {
      RR tmp1, tmp2, tmp3;
      tmp1 = 1;
      tmp2 = val;
      rtn = coeff[0] * tmp1 + coeff[1] * tmp2;
      for (int i = 2; i <= deg; i++) {
        tmp3 = RR(2.0) * val * tmp2 - tmp1;
        tmp1 = tmp2;
        tmp2 = tmp3;
        rtn += coeff[i] * tmp3;
      }
    }
  } else if (type == 2) {
    if (deg == 0) {
      return RR(1);
    } else {
      RR tmp1, tmp2, tmp3;
      RR iden2 = RR(2.0) * val / RR(scale);

      tmp1 = 1;
      tmp2 = val / RR(scale);
      rtn = coeff[0] * tmp1 + coeff[1] * tmp2;
      for (int i = 2; i <= deg; i++) {
        tmp3 = iden2 * tmp2 - tmp1;
        tmp1 = tmp2;
        tmp2 = tmp3;
        rtn += coeff[i] * tmp3;
      }
    }
  }

  return rtn;
}
RR eval(long deg, vector<RR> &coeff, RR val, int type, RR scale) {
  RR rtn;
  if (type == 0) {
    RR tmp;
    tmp = 1;
    rtn = coeff[0] * tmp;
    for (int i = 1; i <= deg; i++) {
      tmp = val * tmp;
      rtn += coeff[i] * tmp;
    }
  } else if (type == 1) {
    if (deg == 0) {
      return RR(1);
    } else {
      RR tmp1, tmp2, tmp3;
      tmp1 = 1;
      tmp2 = val;
      rtn = coeff[0] * tmp1 + coeff[1] * tmp2;
      for (int i = 2; i <= deg; i++) {
        tmp3 = RR(2.0) * val * tmp2 - tmp1;
        tmp1 = tmp2;
        tmp2 = tmp3;
        rtn += coeff[i] * tmp3;
      }
    }
  } else if (type == 2) {
    if (deg == 0) {
      return RR(1);
    } else {
      RR tmp1, tmp2, tmp3;
      RR iden2 = RR(2.0) * val / RR(scale);

      tmp1 = 1;
      tmp2 = val / RR(scale);
      rtn = coeff[0] * tmp1 + coeff[1] * tmp2;
      for (int i = 2; i <= deg; i++) {
        tmp3 = iden2 * tmp2 - tmp1;
        tmp1 = tmp2;
        tmp2 = tmp3;
        rtn += coeff[i] * tmp3;
      }
    }
  }

  return rtn;
}
void showgraph(ofstream &out, RR *coeff, long deg, RR start, RR end, RR sc, int type, RR scale) {
  RR scan;
  scan = start;
  while (scan <= end) {
    out << scan << "," << eval(deg, coeff, scan, type, scale) << endl;
    scan += sc;
  }
}
void showgraph(ofstream &out, vector<RR> coeff, long deg, RR start, RR end, RR sc, int type, RR scale) {
  RR scan;
  scan = start;
  while (scan <= end) {
    out << scan << "," << eval(deg, coeff, scan, type, scale) << endl;
    scan += sc;
  }
}
// deg: odd
// -32 <= expx <= -1.0 && 1.0 < expx <= 32
RR expmaxerr(long deg, RR expx) {
  ifstream E[32];
  long num;
  E[deg].open("result/E" + to_string(deg) + ".txt");
  RR *X, *Y, expy;
  E[deg] >> num;
  X = new RR[num];
  Y = new RR[num];
  for (int i = 0; i < num; i++) {
    E[deg] >> X[i];
    E[deg] >> Y[i];
  }
  for (long i = 1; i < num; i++) {
    if (expx <= X[i]) {
      // -1.0 < expx < 1.0  or  expy ~~ 1.0
      if (X[i - 1] < 0 && X[i] > 0)
        expy = Y[i];
      else if (Y[i - 1] < 0 && Y[i] > 0)
        expy = Y[i];

      // normal case
      else
        expy = Y[i - 1] + (Y[i] - Y[i - 1]) * (expx - X[i - 1]) / (X[i] - X[i - 1]);
      break;
    }
  }
  delete[] X;
  delete[] Y;
  return expy;
}
// deg: odd
// -infty <= expy <= -1.0 && 1.0 < expy < 32
RR invexpmaxerr(long deg, RR expy, RR *X, RR *Y, long num) {
  RR expx;

  for (long i = 1; i < num; i++) {
    if (expy <= Y[i]) {
      // expx ~~ 1.0  or -1.0 < expy < 1.0
      if (X[i - 1] < 0 && X[i] > 0)
        expx = X[i];
      else if (Y[i - 1] < 0 && Y[i] > 0)
        expx = X[i];

      // normal case
      else
        expx = X[i - 1] + (X[i] - X[i - 1]) * (expy - Y[i - 1]) / (Y[i] - Y[i - 1]);
      break;
    }

    // expy > 32.0
    if (i == num - 1) {
      expx = X[i];
      break;
    }
  }
  return expx;
}
RR invexpmaxerr(long deg, RR expy, vector<RR> &X, vector<RR> &Y, long num) {
  RR expx;

  for (long i = 1; i < num; i++) {
    if (expy <= Y[i]) {
      // expx ~~ 1.0  or -1.0 < expy < 1.0
      if (X[i - 1] < 0 && X[i] > 0)
        expx = X[i];
      else if (Y[i - 1] < 0 && Y[i] > 0)
        expx = X[i];

      // normal case
      else
        expx = X[i - 1] + (X[i] - X[i - 1]) * (expy - Y[i - 1]) / (Y[i] - Y[i - 1]);
      break;
    }

    // expy > 32.0
    if (i == num - 1) {
      expx = X[i];
      break;
    }
  }
  return expx;
}
// deg: odd
RR getmaxerr(RR (*func)(RR), vector<RR> &coeff, long deg, RR start, RR end, int type, RR scale, long prec, long num, bool is_opt_sampling) {
  Point *ext = new Point[2 * deg];
  int ext_count;

  // for debugging
  RR sc = (end - start) / static_cast<RR>(num);
  //	find_extreme(sgn, ext, ext_count, coeff, deg, start, end, prec, sc, type, scale);
  find_extreme((*func), ext, ext_count, coeff, deg, start, end, prec, sc, type, scale, is_opt_sampling);
  RR maxerr = RR(0);
  for (long i = 0; i < ext_count; i++) {
    if (maxerr < abs(ext[i].y)) maxerr = abs(ext[i].y);
  }

  return maxerr;
}

// deg: odd
// upgrade version. vector!!!
// input: coeff, deg, start, end, prec, scan, type, scale
// output: ext points, ext_count
RR find_extreme(RR (*func)(RR), Point *&ext, int &ext_count, vector<RR> coeff, long deg, RR start, RR end, long prec, RR scan, int type, RR scale, bool is_opt_sampling) {
  // Parameter setting
  long inc_1 = 0;
  long inc_2 = 0;
  RR scan_1, scan_2, scan_prev;
  RR scan_y1;
  RR scan_y2;
  RR maxerr;
  size_t s = 0;  // this is for optimized sampling
                 //	long scan_prec = 50;
                 //	long temp_prec = scan_prec;

  RR sc = scan;
  RR origin_sc = sc;
  //	RR K = RR(1000);

  long search_inc1, search_inc2, search_inc3, search_inc4, search_iter;
  RR search_start, search_end, search_sc;

  // add the boundary starting point
  ext_count = 0;
  ext[ext_count].x = start;
  ext[ext_count].y = eval(deg, coeff, start, type, scale) - (*func)(start);
  if (ext[ext_count].y > 0)
    ext[ext_count].locmm = 1;
  else
    ext[ext_count].locmm = -1;
  ext_count += 1;

  // start the scan
  if (is_opt_sampling == true) {
    s = 15;
    sc = origin_sc / pow(10, s);
  } else
    sc = origin_sc;

  scan_1 = start;
  scan_2 = start + sc;
  scan_y1 = eval(deg, coeff, scan_1, type, scale) - (*func)(scan_1);
  scan_y2 = eval(deg, coeff, scan_2, type, scale) - (*func)(scan_2);
  if (scan_y1 < scan_y2)
    inc_2 = 1;
  else if (scan_y1 > scan_y2)
    inc_2 = -1;

  while (1) {
    if (scan_2 > end) std::runtime_error("scan2 > end");
    if (is_opt_sampling == true) {
      for (size_t i = 0; i < s; i++) {
        if (start + 10 * origin_sc / pow(10, i) < scan_2 && scan_2 < end - 10 * origin_sc / pow(10, i)) {
          sc = origin_sc / pow(10, i);
          break;
        }
        if (i == s - 1) {
          sc = origin_sc / pow(10, i + 1);
          break;
        }
      }
    } else {
      sc = origin_sc;
    }

    // break condition
    if (scan_2 + sc >= end) break;

    // scan move
    inc_1 = inc_2;
    scan_prev = scan_1;
    scan_1 = scan_2;
    scan_2 = scan_1 + sc;

    scan_y1 = scan_y2;
    //	scan_y2 = eval(deg, coeff, scan_2, type, scale) - sgn(scan_2);
    scan_y2 = eval(deg, coeff, scan_2, type, scale) - (*func)(scan_2);
    if (scan_y1 < scan_y2)
      inc_2 = 1;
    else if (scan_y1 > scan_y2)
      inc_2 = -1;
    else {
      inc_2 = 0;
      for (int i = 0; i < 2; i++) cout << "slope 0 occur!" << endl;
    }

    // binary search
    if (inc_1 == 1 && inc_2 != 1) {
      // search variable initialization
      search_iter = 0;
      search_start = scan_prev;
      search_end = scan_2;
      search_sc = (search_end - search_start) / 4;

      // binary search
      while (search_iter < prec) {
        // obtain the slopes of subintervals
        search_inc1 = (eval(deg, coeff, search_start, type, scale) - (*func)(search_start) < eval(deg, coeff, search_start + search_sc, type, scale) - (*func)(search_start + search_sc) ? 1 : -1);
        search_inc2 = (eval(deg, coeff, search_start + search_sc, type, scale) - (*func)(search_start + search_sc) < eval(deg, coeff, search_start + 2 * search_sc, type, scale) - (*func)(search_start + 2 * search_sc) ? 1 : -1);
        search_inc3 = (eval(deg, coeff, search_start + 2 * search_sc, type, scale) - (*func)(search_start + 2 * search_sc) < eval(deg, coeff, search_end - search_sc, type, scale) - (*func)(search_end - search_sc) ? 1 : -1);
        search_inc4 = (eval(deg, coeff, search_end - search_sc, type, scale) - (*func)(search_end - search_sc) < eval(deg, coeff, search_end, type, scale) - (*func)(search_end) ? 1 : -1);

        // binary search window update
        if (search_inc1 == 1 && search_inc2 == -1) {
          search_end -= 2 * search_sc;
          search_sc /= 2;
        } else if (search_inc2 == 1 && search_inc3 == -1) {
          search_start += search_sc;
          search_end -= search_sc;
          search_sc /= 2;
        } else if (search_inc3 == 1 && search_inc4 == -1) {
          search_start += 2 * search_sc;
          search_sc /= 2;
        }
        search_iter++;
      }

      // add the extreme point
      ext[ext_count].x = (search_start + search_end) / 2;
      //	ext[ext_count].y = eval(deg, coeff, ext[ext_count].x, type, scale) - sgn(ext[ext_count].x);
      ext[ext_count].y = eval(deg, coeff, ext[ext_count].x, type, scale) - (*func)(ext[ext_count].x);
      ext[ext_count].locmm = 1;
      ext_count += 1;
    } else if (inc_1 == -1 && inc_2 != -1) {
      // search variable initialization
      search_iter = 0;
      search_start = scan_prev;
      search_end = scan_2;
      search_sc = (search_end - search_start) / 4;

      // binary search
      while (search_iter < prec) {
        // obtain the slopes of subintervals
        search_inc1 = (eval(deg, coeff, search_start, type, scale) - (*func)(search_start) < eval(deg, coeff, search_start + search_sc, type, scale) - (*func)(search_start + search_sc) ? 1 : -1);
        search_inc2 = (eval(deg, coeff, search_start + search_sc, type, scale) - (*func)(search_start + search_sc) < eval(deg, coeff, search_start + 2 * search_sc, type, scale) - (*func)(search_start + 2 * search_sc) ? 1 : -1);
        search_inc3 = (eval(deg, coeff, search_start + 2 * search_sc, type, scale) - (*func)(search_start + 2 * search_sc) < eval(deg, coeff, search_end - search_sc, type, scale) - (*func)(search_end - search_sc) ? 1 : -1);
        search_inc4 = (eval(deg, coeff, search_end - search_sc, type, scale) - (*func)(search_end - search_sc) < eval(deg, coeff, search_end, type, scale) - (*func)(search_end) ? 1 : -1);

        // binary search window update
        if (search_inc1 == -1 && search_inc2 == 1) {
          search_end -= 2 * search_sc;
          search_sc /= 2;
        } else if (search_inc2 == -1 && search_inc3 == 1) {
          search_start += search_sc;
          search_end -= search_sc;
          search_sc /= 2;
        } else if (search_inc3 == -1 && search_inc4 == 1) {
          search_start += 2 * search_sc;
          search_sc /= 2;
        }
        search_iter++;
      }

      // add the extreme point
      ext[ext_count].x = (search_start + search_end) / 2;
      //	ext[ext_count].y = eval(deg, coeff, ext[ext_count].x, type, scale) - sgn(ext[ext_count].x);
      ext[ext_count].y = eval(deg, coeff, ext[ext_count].x, type, scale) - (*func)(ext[ext_count].x);
      ext[ext_count].locmm = -1;
      ext_count += 1;
    }
  }

  // add the ending boundary point
  ext[ext_count].x = end;
  //	ext[ext_count].y = eval(deg, coeff, end, type, scale) - sgn(end);
  ext[ext_count].y = eval(deg, coeff, end, type, scale) - (*func)(end);
  if (ext[ext_count].y > 0)
    ext[ext_count].locmm = 1;
  else
    ext[ext_count].locmm = -1;
  ext_count += 1;

  // calculate maximum error maxerr
  maxerr = 0;
  for (long i = 0; i < ext_count; i++) {
    if (maxerr < abs(ext[i].y)) maxerr = abs(ext[i].y);
  }

  return maxerr;
}
}  // namespace minicomp
