#pragma once
#include <iostream>
#include <vector>

#include "ckks_evaluator.cuh"
#include "phantom.h"

using namespace std;
using namespace phantom;
using namespace nexus;

namespace nexus {
class ArgmaxEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  ArgmaxEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}

  void argmax(PhantomCiphertext &x, PhantomCiphertext &res, int len);
};
}  // namespace nexus
