#pragma once

#include <seal/seal.h>

#include <vector>

#include "ckks_evaluator.h"

using namespace std;
using namespace seal;

class GeLUEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  GeLUEvaluator(CKKSEvaluator &ckks) {
    this->ckks = &ckks;
  }

  void gelu(Ciphertext &x, Ciphertext &res);
  vector<double> gelu_plain(vector<double> &input);
};
