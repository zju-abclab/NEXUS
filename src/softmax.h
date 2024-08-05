#include <seal/seal.h>

#include <iostream>
#include <vector>

#include "ckks_evaluator.h"

class SoftmaxEvaluator {
 public:
  CKKSEvaluator *ckks = nullptr;
  SoftmaxEvaluator(CKKSEvaluator &ckks) {
    this->ckks = &ckks;
  }
  void softmax(Ciphertext &x, Ciphertext &res, int len);
};
