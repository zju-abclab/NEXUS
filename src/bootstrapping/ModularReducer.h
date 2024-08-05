#pragma once

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

#include "Polynomial.h"
#include "RemezArcsin.h"
#include "RemezCos.h"
// #include "ScaleInvEvaluator.h"

using namespace std;

class ModularReducer {
 public:
  long boundary_K;
  double log_width;
  long deg;
  long num_double_formula;

  double inverse_log_width;
  long inverse_deg;

  double scale_inverse_coeff;

  SEALContext &context;
  CKKSEncoder &encoder;
  Encryptor &encryptor;
  Evaluator &evaluator;
  RelinKeys &relin_keys;
  Decryptor &decryptor;

  RemezParam rmparm;

  RemezCos *poly_generator;
  RemezArcsin *inverse_poly_generator;

  boot::Polynomial sin_cos_polynomial;
  boot::Polynomial inverse_sin_polynomial;

  ModularReducer(long _boundary_K, double _log_width, long _deg, long _num_double_formula, long _inverse_deg,
                 SEALContext &_context,
                 CKKSEncoder &_encoder,
                 Encryptor &_encryptor,
                 Evaluator &_evaluator,
                 RelinKeys &_relin_keys,
                 Decryptor &_decryptor);

  void double_angle_formula(Ciphertext &cipher);
  void double_angle_formula_scaled(Ciphertext &cipher, double scale_coeff);
  void generate_sin_cos_polynomial();
  void generate_inverse_sine_polynomial();
  void write_polynomials();
  void modular_reduction(Ciphertext &rtn, Ciphertext &cipher);
};
