#pragma once

#include <NTL/RR.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "func.h"
// #include"ScaleInvEvaluator.h"

using namespace std;
using namespace NTL;
using namespace seal;
using namespace seal::util;

namespace boot {
class Polynomial {
 public:
  RR *coeff = 0;
  long deg, heap_k, heap_m, heaplen;
  RR *chebcoeff = 0;
  Polynomial **poly_heap = 0;

  Polynomial();
  Polynomial(long _deg);
  Polynomial(long _deg, RR *_coeff, string tag);

  ~Polynomial();
  void set_polynomial(long _deg, RR *_coeff, string tag);
  void set_zero_polynomial(long _deg);
  void showcoeff();
  void showchebcoeff();
  void copy(Polynomial &poly);
  void power_to_cheb();
  void cheb_to_power();
  RR evaluate(RR value);

  void constmul(RR constant);
  void mulinplace(Polynomial &poly);
  void addinplace(Polynomial &poly);
  void subtinplace(Polynomial &poly);

  void change_variable_scale(RR scale);  // f_new(x) = f(x / scale)
  void generate_poly_heap_manual(long k, long m);
  void generate_poly_heap();
  void generate_poly_heap_odd();

  void write_heap_to_file(ofstream &out);
  void read_heap_from_file(ifstream &in);

  // void homomorphic_poly_evaluation(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, ScaleInvEvaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor);
  void homomorphic_poly_evaluation(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, Evaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor);
};

void mul(Polynomial &rtn, Polynomial &a, Polynomial &b);
void add(Polynomial &rtn, Polynomial &a, Polynomial &b);
void subt(Polynomial &rtn, Polynomial &a, Polynomial &b);

void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider);
void chebyshev(Polynomial &rtn, long deg);
void second_chebyshev_times_x_for_sine(Polynomial &rtn, long deg);
}  // namespace boot
