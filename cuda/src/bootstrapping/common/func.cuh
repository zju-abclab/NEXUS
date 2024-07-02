#pragma once

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#include "Point.cuh"

using namespace std;
using namespace NTL;

RR fracpart(RR x);
RR sqrtfracpart(RR x, RR a);
RR fraccos(RR x, long scale);
RR arcsin(RR x);
bool yabscompare(Point a, Point b);
bool ycompare(Point a, Point b);
bool xcompare(Point a, Point b);
bool isin(RR x, long K, RR width);
RR chebeval(long deg, RR *coeff, RR val);
void showgraph(ofstream &out, RR *coeff, long deg, long K, RR sc);
bool oddevennextcombi(long *arr, long arrlen, long len);
void oddbabycount(long &k, long &m, long deg);
void babycount(long &k, long &m, long deg);

void add(complex<double> *&rtn, complex<double> *&vec1, complex<double> *&vec2, long n);
void addinplace(complex<double> *&vec, complex<double> *&addvec, long n);
void subt(complex<double> *&rtn, complex<double> *&vec1, complex<double> *&vec2, long n);
void subtinplace(complex<double> *&vec, complex<double> *&subtvec, long n);
void mul(complex<double> *&rtn, complex<double> *&vec1, complex<double> *&vec2, long n);
void mulinplace(complex<double> *&vec, complex<double> *&mulvec, long n);
void constmul(complex<double> *&rtn, complex<double> *&vec, complex<double> constant, long n);
void constmulinplace(complex<double> *&vec, complex<double> constant, long n);
void text_to_array(ifstream &in, RR *&array, long n);
int giantstep(int M);
void rotation(int logslot, int Nh, int shiftcount, const vector<complex<double>> &vec, vector<complex<double>> &rtnvec);

int max_index(double *array, int length);
