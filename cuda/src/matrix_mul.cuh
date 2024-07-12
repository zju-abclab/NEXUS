#pragma once
#include <fstream>
#include <iostream>
#include <vector>

#include "ckks_evaluator.cuh"
#include "phantom.h"
#include "troy.cuh"

namespace nexus {
using namespace std;
using namespace phantom;
using namespace troy;

class MMEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;
  int print_count = 0;

  void enc_compress_ciphertext(vector<double> &values, PhantomCiphertext &ct);
  void enc_compress_ciphertext(vector<double> &values, Ciphertext &ct);

  vector<PhantomCiphertext> decompress_ciphertext(const PhantomCiphertext &encrypted);
  vector<Ciphertext> decompress_ciphertext(const Ciphertext &encrypted);

 public:
  MMEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}

  // Helper functions
  inline vector<vector<double>> read_matrix(const std::string &filename, int rows, int cols) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    std::ifstream file(filename);

    if (!file.is_open()) {
      std::cerr << "Can not open file: " << filename << std::endl;
      return matrix;
    }

    std::string line;
    for (int i = 0; i < rows; ++i) {
      if (std::getline(file, line)) {
        std::istringstream iss(line);
        for (int j = 0; j < cols; ++j) {
          if (!(iss >> matrix[i][j])) {
            std::cerr << "read error: " << filename << " (row: " << i << ", column: " << j << ")" << std::endl;
          }
        }
      }
    }

    file.close();
    return matrix;
  }

  inline vector<vector<double>> transpose_matrix(const vector<vector<double>> &matrix) {
    if (matrix.empty()) {
      return {};
    }
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<double>> transposedMatrix(cols, std::vector<double>(rows));

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        transposedMatrix[j][i] = matrix[i][j];
      }
    }

    return transposedMatrix;
  }

  // Evaluation function
  void matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<PhantomCiphertext> &res);
  void matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<Ciphertext> &res);

  void multiply_power_of_x(PhantomCiphertext &encrypted, PhantomCiphertext &destination, int index);
  void multiply_power_of_x(Ciphertext &encrypted, Ciphertext &destination, int index);
};
}  // namespace nexus
