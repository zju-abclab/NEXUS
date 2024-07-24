#include <chrono>
#include <iostream>
#include <vector>

#include "Bootstrapper.cuh"
#include "ckks_evaluator.cuh"
#include "phantom.h"

using namespace std;
using namespace chrono;
using namespace phantom;

class ArgmaxEvaluator {
 private:
  // Bootstrapping
  long boundary_K = 25;
  long deg = 59;
  long scale_factor = 2;
  long inverse_deg = 1;
  long loge = 10;

  int bs_mod_count = 14;

  long logN;
  long logn;
  int total_level;
  long sparse_slots;

 public:
  CKKSEvaluator *ckks = nullptr;

  // execution time minus bootstrapping setup time
  chrono::duration<int64_t, nano> time_elapsed = chrono::duration<int64_t, nano>(0);

  ArgmaxEvaluator(CKKSEvaluator &ckks, int main_mod_count) {
    this->ckks = &ckks;

    this->logN = log2(ckks.degree);
    this->logn = logN - 2;
    this->sparse_slots = (1 << logn);

    this->total_level = main_mod_count + bs_mod_count;
  }

  void argmax(PhantomCiphertext &x, PhantomCiphertext &res, int len);

  void bootstrap(PhantomCiphertext &x);
};
