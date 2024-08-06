#pragma once

#include <seal/seal.h>

#include "Bootstrapper.h"
#include "ckks_evaluator.h"

using namespace std;
using namespace seal;

class ArgmaxEvaluator {
 private:
  Bootstrapper *bootstrapper = nullptr;

 public:
  CKKSEvaluator *ckks = nullptr;

  ArgmaxEvaluator(CKKSEvaluator &ckks, Bootstrapper &bootstrapper) {
    this->ckks = &ckks;
    this->bootstrapper = &bootstrapper;
  }

  void argmax(Ciphertext &x, Ciphertext &res, int len);

  void bootstrap(Ciphertext &x);
};
