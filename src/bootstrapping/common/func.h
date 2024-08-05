#pragma once

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#include "Point.h"
#include "seal/seal.h"

using namespace std;
using namespace NTL;
using namespace seal;

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

void make_modulus_equal(shared_ptr<SEALContext> &context, Evaluator &evaluator, double scale, Ciphertext &cipher1, Ciphertext &cipher2, Ciphertext &rtncipher1, Ciphertext &rtncipher2);
int max_index(double *array, int length);
void decrypt_and_print(const Ciphertext &cipher, Decryptor &decryptor, CKKSEncoder &encoder, long sparse_slots, size_t front = 5, size_t back = 5);
void decrypt_and_print_and_max_round(const Ciphertext &cipher, Decryptor &decryptor, CKKSEncoder &encoder, double unit, long sparse_slots, size_t front = 5, size_t back = 5);
void decode_and_print(const Plaintext &plain, CKKSEncoder &encoder, long sparse_slots, size_t front = 5, size_t back = 5);
void decrypt_and_check_mod(const Ciphertext &mod_input, const Ciphertext &mod_result, Decryptor &decryptor, CKKSEncoder &encoder, double unit, long sparse_slots, size_t front = 5, size_t back = 5);
