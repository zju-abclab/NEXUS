// #include "Remez.h"
#include "MinicompRemez.cuh"

namespace minicomp {
Remez::Remez(RR (*_func)(RR), size_t _inter_num, vector<RR> _inter_start, vector<RR> _inter_end, RR _sc, long _prec, long _deg, long _iter, int _type, RR _scale, bool _is_opt_sampling) {
  sc = _sc;
  prec = _prec;
  deg = _deg;
  iter = _iter;
  type = _type;
  scale = _scale;
  inter_num = _inter_num;
  inter_start = _inter_start;
  inter_end = _inter_end;
  func = _func;
  is_opt_sampling = _is_opt_sampling;

  sam = new Point[deg + 2];
  coeff = vector<RR>(deg + 1, RR(0));
  ext = new Point[2 * deg];
  ext_maxsum_index = new int[deg + 2];
  temp_ext = new RR[2 * deg];

  v.SetLength(deg + 2);
  w.SetLength(deg + 2);
  v_0.SetLength(deg + 2);
  m.SetDims(deg + 2, deg + 2);
}

Remez::~Remez() {
  //	cout << "~Remez" << endl;
  delete[] sam;
  //	delete [] coeff;
  delete[] ext;
  delete[] ext_maxsum_index;
  delete[] temp_ext;
  //	cout << "~Remez finish" << endl;
}

void Remez::initialize() {
  size_t n, m, remain, t;
  n = deg + 1;
  t = inter_num;
  m = (n + 1 + t - 1) / t;  // ceiling of (n+1)/t
  remain = n + 1 - m * (t - 1);
  for (size_t i = 0; i < t - 1; i++) {
    for (size_t j = 0; j < m; j++) {
      sam[m * i + j].x = inter_start[i] + (inter_end[i] - inter_start[i]) / static_cast<RR>(m + 1) * static_cast<RR>(j + 1);
      sam[m * i + j].y = (*func)(sam[m * i + j].x);
    }
  }
  for (size_t j = 0; j < remain; j++) {
    sam[m * (t - 1) + j].x = inter_start[t - 1] + (inter_end[t - 1] - inter_start[t - 1]) / static_cast<RR>(remain + 1) * static_cast<RR>(j + 1);
    sam[m * (t - 1) + j].y = (*func)(sam[m * (t - 1) + j].x);
  }

  // cout << "initialize points" << endl;
  // for(int i=0; i<n+1; i++) cout << sam[i].x << ", " << sam[i].y << endl;
}
void Remez::getcoeffwitherr() {
  for (long i = 0; i < deg + 2; i++) {
    v[i] = sam[i].x;
    w[i] = sam[i].y;
  }
  RR var;

  // powereval
  if (type == 0) {
    for (long i = 0; i < deg + 2; i++) {
      var = v[i];
      m[i][0] = 1;
      for (long j = 1; j < deg + 1; j++) {
        m[i][j] = var * m[i][j - 1];
      }
      m[i][deg + 1] = 2 * (i % 2) - 1;
    }
  }
  // chebeval
  else if (type == 1) {
    for (long i = 0; i < deg + 2; i++) {
      var = v[i];
      m[i][0] = 1;
      m[i][1] = var;
      for (long j = 2; j < deg + 1; j++) {
        m[i][j] = RR(2.0) * var * m[i][j - 1] - m[i][j - 2];
      }
      m[i][deg + 1] = 2 * (i % 2) - 1;
    }
  } else if (type == 2) {
    for (long i = 0; i < deg + 2; i++) {
      var = v[i];
      m[i][0] = 1;
      m[i][1] = var / RR(scale);
      for (long j = 2; j < deg + 1; j++) {
        m[i][j] = RR(2.0) * var / RR(scale) * m[i][j - 1] - m[i][j - 2];
      }
      m[i][deg + 1] = 2 * (i % 2) - 1;
    }
  }

  RR deter;
  mat_RR mtrns;
  mat_RR minv;

  transpose(mtrns, m);
  inv(deter, minv, mtrns);
  v_0 = w * minv;

  for (int i = 0; i <= deg; i++) {
    coeff[i] = v_0[i];
  }
}

int Remez::getextreme() {
  // find extreme points
  Point* temp_ext = new Point[2 * deg];
  int temp_ext_count;
  ext_count = 0;

  for (size_t j = 0; j < inter_num; j++) {
    find_extreme(func, temp_ext, temp_ext_count, coeff, deg, inter_start[j], inter_end[j], prec, sc, type, scale, is_opt_sampling);
    for (int i = 0; i < temp_ext_count; i++) {
      ext[ext_count].x = temp_ext[i].x;
      ext[ext_count].y = temp_ext[i].y;
      ext[ext_count].locmm = temp_ext[i].locmm;
      ext_count++;
    }
  }
  // cout << "number of extreme points: " << ext_count << endl;
  // cout << "print extreme points" << endl;
  // for(size_t i=0; i<ext_count; i++) cout << ext[i].x << ", " << ext[i].y << endl;

  /*
  //	find_extreme(temp_ext, temp_ext_count, coeff, deg, first_inter_start, first_inter_end, prec, sc, type, scale);
          find_extreme(temp_ext, temp_ext_count, coeff, deg, inter_start[0], inter_end[0], prec, sc, type, scale);
          for(int i=0; i<temp_ext_count; i++) {
                  ext[ext_count].x = temp_ext[i].x;
                  ext[ext_count].y = temp_ext[i].y;
                  ext[ext_count].locmm = temp_ext[i].locmm;
                  ext_count ++;
          }

  //	find_extreme(temp_ext, temp_ext_count, coeff, deg, second_inter_start, second_inter_end, prec, sc, type, scale);
          find_extreme(temp_ext, temp_ext_count, coeff, deg, inter_start[1], inter_end[1], prec, sc, type, scale);
          for(int i=0; i<temp_ext_count; i++) {
                  ext[ext_count].x = temp_ext[i].x;
                  ext[ext_count].y = temp_ext[i].y;
                  ext[ext_count].locmm = temp_ext[i].locmm;
                  ext_count ++;
          }
  */

  // show error message
  if (ext_count < deg + 2) {
    cout << "Couldn't find all the extreme points! " << endl;
    cout << "ext_count: " << ext_count << endl;
    return -1;
  }
  return 1;
}

void Remez::choosemaxs() {
  for (int i = 0; i < ext_count; i++) temp_ext[i] = abs(ext[i].y);
  MaxSubsetSum(temp_ext, ext_count, deg + 2, ext_maxsum_index);

  maxerr = RR(0);
  for (int i = 0; i < deg + 2; i++) {
    sam[i].x = ext[ext_maxsum_index[i]].x;
    sam[i].y = (*func)(ext[ext_maxsum_index[i]].x);
    if (maxerr < ext[ext_maxsum_index[i]].y) maxerr = ext[ext_maxsum_index[i]].y;
  }
}

bool Remez::test() {
  long it = 0;
  // int sum = 0;
  // RR* temp = new RR[deg+2];

  initialize();
  while (it < iter) {
    it++;
    cout << "iteration: " << it << endl;
    getcoeffwitherr();
    if (getextreme() < 0) return false;
    choosemaxs();

    RR max, min;

    max = abs(ext[0].y);
    min = abs(ext[0].y);
    for (int i = 1; i < ext_count; i++) {
      if (abs(ext[i].y) > max) max = abs(ext[i].y);
      if (abs(ext[i].y) < min) min = abs(ext[i].y);
    }
    if ((max - min) / min < pow(2.0, -60.0)) {
      cout << "delta: " << (max - min) / min << endl;
      break;
    }
  }

  // delete [] temp;
  return true;
}

void Remez::printgraph(ofstream& out, RR start, RR end, RR sc) {
  RR scan;
  scan = start;
  while (scan <= end) {
    out << scan << " " << eval(deg, coeff, scan, type, scale) << endl;
    scan += sc;
  }
}

RR Remez::getmaxerr() {
  return maxerr;
}

void Remez::getcoeff(RR _coeff[]) {
  for (long i = 0; i < deg + 1; i++)
    _coeff[i] = coeff[i];
}

void Remez::getcoeff(vector<RR>& _coeff) {
  _coeff.clear();
  for (long i = 0; i < deg + 1; i++)
    _coeff.emplace_back(coeff[i]);
}

void Remez::getext_xpos(vector<RR>& _ext_xpos) {
  _ext_xpos.clear();
  for (int i = deg / 2 + 1; i < deg + 2; i++) _ext_xpos.emplace_back(sam[i].x);
}
}  // namespace minicomp
