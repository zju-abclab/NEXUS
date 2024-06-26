#pragma once
#include "ckks_evaluator.cuh"
#include "phantom.h"

using namespace std;
using namespace phantom;
using namespace nexus;

namespace nexus {
class LNEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  LNEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}
  void layer_norm(PhantomCiphertext &x, PhantomCiphertext &res, int len);
};
}  // namespace nexus
