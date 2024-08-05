#include <seal/ciphertext.h>

#include <vector>

#include "ckks_evaluator.h"
class MMEvaluator {
 private:
  vector<Ciphertext> expand_ciphertext(const Ciphertext &encrypted, uint32_t m, GaloisKeys &galkey, vector<uint32_t> &galois_elts);

  void multiply_power_of_x(Ciphertext &encrypted, Ciphertext &destination, int index);

  void enc_compress_ciphertext(vector<double> &vec, Ciphertext &ct);

 public:
  CKKSEvaluator *ckks = nullptr;
  MMEvaluator(CKKSEvaluator &ckks) {
    this->ckks = &ckks;
  }

  void matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<Ciphertext> &res);

  std::vector<std::vector<double>> readMatrix(const std::string &filename, int rows, int cols);

  std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>> &matrix);
};
