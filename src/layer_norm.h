#pragma once

#include <seal/seal.h>

#include "ckks_evaluator.h"

using namespace std;
using namespace seal;

class LNEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  LNEvaluator(CKKSEvaluator &ckks) {
    this->ckks = &ckks;
  }

  void layer_norm(Ciphertext &x, Ciphertext &res, int len);
};
