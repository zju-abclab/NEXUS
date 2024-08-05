#include <seal/seal.h>

#include <iostream>
#include <vector>

#include "ckks_evaluator.h"

class LNEvaluator {
 public:
  CKKSEvaluator *ckks = nullptr;
  LNEvaluator(CKKSEvaluator &ckks) {
    this->ckks = &ckks;
  }
  void layer_norm(Ciphertext &x, Ciphertext &res, int len);
};
