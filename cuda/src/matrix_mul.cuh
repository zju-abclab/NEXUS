#pragma once
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include "ckks_evaluator.cuh"
#include "kernels.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;
using namespace phantom::util;

class MMEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

  void expand_encode(vector<double> &vec, PhantomCiphertext &ct);

  // void multiply_power_of_x(PhantomCiphertext &encrypted, PhantomCiphertext &destination, int index);
  vector<PhantomCiphertext> expand_ciphertext(const PhantomCiphertext &encrypted, uint32_t m, PhantomGaloisKey &galkey, vector<uint32_t> &galois_elts);

  void expand_and_insert(CKKSEvaluator *ckks, mutex &mtx, const PhantomCiphertext &compressed_ct, vector<PhantomCiphertext> &b_expanded_cts);
  void multithreaded_expansion(CKKSEvaluator *ckks, const vector<PhantomCiphertext> &b_compressed_cts, vector<PhantomCiphertext> &b_expanded_cts);

 public:
  MMEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}

  // Helper functions
  vector<vector<double>> read_matrix(const std::string &filename, int rows, int cols);
  vector<vector<double>> transpose_matrix(const vector<vector<double>> &matrix);

  // Evaluation function
  void matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<PhantomCiphertext> &res);
  void multiply_power_of_x(PhantomCiphertext &encrypted, PhantomCiphertext &destination, int index);
};
}  // namespace nexus
