#pragma once

#include <seal/seal.h>

#include "ckks_evaluator.h"

using namespace std;
using namespace seal;

class SoftmaxEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  SoftmaxEvaluator(CKKSEvaluator &ckks) {
    this->ckks = &ckks;
  }

  void softmax(Ciphertext &x, Ciphertext &res, int len);
};
