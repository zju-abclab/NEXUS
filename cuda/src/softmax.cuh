#pragma once

#include "ckks_evaluator.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;

class SoftmaxEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  SoftmaxEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}

  void softmax(PhantomCiphertext &x, PhantomCiphertext &res, int len);
};
}  // namespace nexus
