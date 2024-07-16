#pragma once

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

#include "RemezArcsin.cuh"
#include "RemezCos.cuh"
#include "common/Polynomial.cuh"

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

  CKKSEvaluator *ckks = nullptr;

  RemezParam rmparm;

  RemezCos *poly_generator;
  RemezArcsin *inverse_poly_generator;

  boot::Polynomial sin_cos_polynomial;
  boot::Polynomial inverse_sin_polynomial;

  ModularReducer(long _boundary_K, double _log_width, long _deg, long _num_double_formula, long _inverse_deg, CKKSEvaluator *ckks);

  void double_angle_formula(PhantomCiphertext &cipher);
  void double_angle_formula_scaled(PhantomCiphertext &cipher, double scale_coeff);
  void generate_sin_cos_polynomial();
  void generate_inverse_sine_polynomial();
  void write_polynomials();
  void modular_reduction(PhantomCiphertext &rtn, PhantomCiphertext &cipher);
};
