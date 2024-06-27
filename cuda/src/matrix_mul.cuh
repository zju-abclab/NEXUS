#pragma once
#include <iostream>
#include <vector>

#include "ckks_evaluator.cuh"
#include "phantom.h"

using namespace std;
using namespace phantom;
using namespace nexus;

namespace nexus {
class MMEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

  vector<PhantomCiphertext> expand_ciphertext(const PhantomCiphertext &encrypted, uint32_t m, PhantomGaloisKey &galkey, vector<uint32_t> &galois_elts);

  void multiply_power_of_x(PhantomCiphertext &encrypted, PhantomCiphertext &destination, int index);

  void expand_encode(vector<double> &vec, PhantomCiphertext &ct);

 public:
  MMEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}

  void matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<PhantomCiphertext> &res);
};
}  // namespace nexus
