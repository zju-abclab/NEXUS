#include "func.cuh"

RR fracpart(RR x) {
  return x - round(x);
}

RR sqrtfracpart(RR x, RR a) {
  ZZ roundx;
  roundx = RoundToZZ(x);
  if (roundx % 2 == 0)
    return sqrt(x - round(x) + a);
  else
    return (-1) * sqrt(x - round(x) + a);
}

RR fraccos(RR x, long scale) {
  RR sc_fac = RR(1 << scale);
  return cos(2 * ComputePi_RR() * (x - 0.25) / sc_fac);
}

RR arcsin(RR x) {
  RR rtn = x;
  while (abs(sin(rtn) - x) > pow(2, -300)) {
    rtn = rtn - (sin(rtn) - x) / cos(rtn);
  }

  return rtn;
}

bool yabscompare(Point a, Point b) {
  return (abs(a.y) > abs(b.y));
}

bool ycompare(Point a, Point b) {
  return (a.y > b.y);
}

bool xcompare(Point a, Point b) {
  return (a.x < b.x);
}

bool isin(RR x, long K, RR width) {
  if (x < -K || x > K)
    return false;
  else if (abs(fracpart(x)) <= width)
    return true;
  else
    return false;
}

RR chebeval(long deg, RR *coeff, RR val) {
  RR tmp1, tmp2, tmp3;
  RR rtn;
  tmp1 = 1;
  tmp2 = val;
  rtn = coeff[0] * tmp1 + coeff[1] * tmp2;
  for (int i = 2; i <= deg; i++) {
    tmp3 = 2 * val * tmp2 - tmp1;
    tmp1 = tmp2;
    tmp2 = tmp3;
    rtn += coeff[i] * tmp3;
  }
  return rtn;
}

void showgraph(ofstream &out, RR *coeff, long deg, long K, RR sc) {
  RR scan;
  scan = (-1) * K;
  while (scan <= K) {
    out << scan << "," << chebeval(deg, coeff, scan / K) << endl;
    scan += sc;
  }
}

bool oddevennextcombi(long *arr, long arrlen, long len) {
  long i = arrlen - 1;
  long ind = len - 1;
  while (ind - arr[i] < 2) {
    i--;
    ind--;
  }

  if (i >= 0) {
    arr[i] += 2;
    i++;
    while (i < arrlen) {
      arr[i] = arr[i - 1] + 1;
      i++;
    }
    return false;
  } else
    return true;
}

void oddbabycount(long &mink, long &minm, long deg) {
  long m, mind = 0;
  long mineval = 100000;
  for (long k = 2; k <= deg; k += 2) {
    m = 1;
    while ((1 << m) * k < deg) {
      m++;
    }
    if (ceil(log(k) / log(2)) + m == ceil(log(deg) / log(2))) {
      if (mineval > ceil((deg + 0.0) / (k + 0.0)) + k / 2 - 5 + m + ceil(log(k)) + max(3.0 - m, ceil((deg + 0.0) / (1.5 * (1 << (m - 1)) * k)))) {
        mineval = ceil((deg + 0.0) / (k + 0.0)) + k / 2 - 5 + m + ceil(log(k)) + max(3.0 - m, ceil((deg + 0.0) / (1.5 * (1 << (m - 1)) * k)));
        mink = k;
        minm = m;
        mind = ceil(log(k) / log(2)) + m;
      } else if (mineval == ceil((deg + 0.0) / (k + 0.0)) + k / 2 - 5 + m + ceil(log(k)) + max(3.0 - m, ceil((deg + 0.0) / (1.5 * (1 << (m - 1)) * k)))) {
        if (mind > ceil(log(k) / log(2)) + m) {
          mink = k;
          minm = m;
          mind = ceil(log(k) / log(2)) + m;
        }
      }
    }
  }
}

void babycount(long &mink, long &minm, long deg) {
  int curr_mul, min_mul;
  mink = 2;

  double d_over_k = static_cast<double>(deg) / mink;
  int log2_d_over_k = ceil(log2(d_over_k));
  int ceil_d_over_k = ceil(d_over_k);
  min_mul = log2_d_over_k + mink + ceil_d_over_k - 3;
  minm = log2_d_over_k;

  for (int i = 3; i < 2 * sqrt(deg); i++) {
    d_over_k = static_cast<double>(deg) / i;
    log2_d_over_k = ceil(log2(d_over_k));
    ceil_d_over_k = ceil(d_over_k);
    curr_mul = log2_d_over_k + i + ceil_d_over_k - 3;

    if (min_mul > curr_mul) {
      mink = i;
      min_mul = curr_mul;
      minm = log2_d_over_k;
    }
  }
}

void add(complex<double> *&rtn, complex<double> *&vec1, complex<double> *&vec2, long n) {
  for (int i = 0; i < n; i++) {
    rtn[i] = vec1[i] + vec2[i];
  }
}

void addinplace(complex<double> *&vec, complex<double> *&addvec, long n) {
  for (int i = 0; i < n; i++) {
    vec[i] += addvec[i];
  }
}
void subt(complex<double> *&rtn, complex<double> *&vec1, complex<double> *&vec2, long n) {
  for (int i = 0; i < n; i++) {
    rtn[i] = vec1[i] - vec2[i];
  }
}

void subtinplace(complex<double> *&vec, complex<double> *&subtvec, long n) {
  for (int i = 0; i < n; i++) {
    vec[i] -= subtvec[i];
  }
}

void mul(complex<double> *&rtn, complex<double> *&vec1, complex<double> *&vec2, long n) {
  for (int i = 0; i < n; i++) {
    rtn[i] = vec1[i] * vec2[i];
  }
}

void mulinplace(complex<double> *&vec, complex<double> *&mulvec, long n) {
  for (int i = 0; i < n; i++) {
    vec[i] *= mulvec[i];
  }
}

void constmul(complex<double> *&rtn, complex<double> *&vec, complex<double> constant, long n) {
  for (int i = 0; i < n; i++) {
    rtn[i] = constant * vec[i];
  }
}

void constmulinplace(complex<double> *&vec, complex<double> constant, long n) {
  for (int i = 0; i < n; i++) {
    vec[i] *= constant;
  }
}

void text_to_array(ifstream &in, RR *&array, long n) {
  long i = 0;
  if (in.is_open()) {
    while (!in.eof() && i < n) {
      std::cout << i << std::endl;
      in >> array[i];
      std::cout << array[i] << std::endl;
      i++;
    }
    std::cout << "ok" << std::endl;
  }
}

int giantstep(int M) {
  int minval = M, mink = 1, currval;
  for (int k = 1; k <= 3 * sqrt(M); k++) {
    currval = ceil((M + 0.0) / (k + 0.0)) + k - 1;
    if (currval < minval) {
      minval = currval;
      mink = k;
    }
  }
  return mink;
}

void rotation(int logslot, int Nh, int shiftcount, const vector<complex<double>> &vec, vector<complex<double>> &rtnvec) {
  int slotlen = (1 << logslot);
  int repeatcount = Nh / slotlen;
  rtnvec.clear();
  for (int j = 0; j < repeatcount; j++) {
    for (int i = 0; i < slotlen; i++) {
      rtnvec.push_back(vec[(slotlen + i + shiftcount) % slotlen]);
    }
  }
}

int max_index(double *array, int length) {
  int max_ind = 0;
  double max = array[0];

  for (int i = 1; i < length; i++) {
    if (array[i] > max) {
      max_ind = i;
      max = array[i];
    }
  }
  return max_ind;
}
