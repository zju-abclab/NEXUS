#include <chrono>
#include <iostream>
#include <vector>

#include "bootstrapping/Bootstrapper.cuh"
#include "ckks_evaluator.cuh"
#include "phantom.h"

using namespace std;
using namespace chrono;
using namespace phantom;

class ArgmaxEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;
  Bootstrapper *bootstrapper = nullptr;

  // execution time minus bootstrapping setup time
  chrono::duration<int64_t, nano> time_elapsed = chrono::duration<int64_t, nano>(0);

 public:
  ArgmaxEvaluator(CKKSEvaluator &ckks, Bootstrapper &bootstrapper, int main_mod_count) {
    this->ckks = &ckks;
    this->bootstrapper = &bootstrapper;
  }

  void argmax(PhantomCiphertext &x, PhantomCiphertext &res, int len);

  void bootstrap(PhantomCiphertext &x);
};
