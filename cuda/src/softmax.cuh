#pragma once
#include <iostream>
#include <vector>

#include "ckks_evaluator.cuh"
#include "phantom.h"

using namespace std;
using namespace phantom;
using namespace nexus;

namespace nexus {
class SoftmaxEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  SoftmaxEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}

  void softmax2(PhantomCiphertext &x, PhantomCiphertext &res, int len);
};
}  // namespace nexus
