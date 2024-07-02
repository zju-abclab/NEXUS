#pragma once

#include <NTL/RR.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "Point.cuh"

using namespace std;

namespace minicomp {
ZZ mod(ZZ a, ZZ q);
long pmod(long a, long b);
long pow2(long n);
size_t ceil_divide(size_t x, size_t y);
size_t ceil_to_int(double x);
int floor_to_int(double x);
long log2_long(long n);
long num_one(long n);
RR GetApproxError(int d, RR t, bool is_first_function);
RR GetInvApproxError(int d, RR t, RR *X, RR *Y, long num);
RR GetInvApproxError(int d, RR t, vector<RR> &X, vector<RR> &Y, long num);
RR real(RR p);
int dep(list<int> &lt);
int dep(int deg);
int mult(int deg);
int enough_mult(int a);
int enough_dep(int a);
RR exptoreal(RR x);
RR realtoexp(RR x);
RR sgn(RR x);
RR comp(RR a, RR b);
RR fracpart(RR x);
RR ReLU(RR x);
double comp(double a, double b);
double ReLU(double x);
bool yabscompare(Point a, Point b);
bool ycompare(Point a, Point b);
bool xcompare(Point a, Point b);
bool isin(RR x, long K, RR width);
RR eval(long deg, RR *coeff, RR val, int type, RR scale);
RR eval(long deg, vector<RR> &coeff, RR val, int type, RR scale);
void showgraph(ofstream &out, RR *coeff, long deg, RR start, RR end, RR sc, int type, RR scale);
void showgraph(ofstream &out, vector<RR> coeff, long deg, RR start, RR end, RR sc, int type, RR scale);
bool oddevennextcombi(long *arr, long arrlen, long len);
RR expmaxerr(long deg, RR expx);
RR invexpmaxerr(long deg, RR expy, RR *X, RR *Y, long num);
RR invexpmaxerr(long deg, RR expy, vector<RR> &X, vector<RR> &Y, long num);
RR getmaxerr(RR (*func)(RR), vector<RR> &coeff, long deg, RR start, RR end, int type, RR scale, long prec, long num, bool is_opt_sampling);
RR find_extreme(RR (*func)(RR), Point *&ext, int &ext_count, vector<RR> coeff, long deg, RR start, RR end, long prec, RR scan, int type, RR scale, bool is_opt_sampling);
}  // namespace minicomp
