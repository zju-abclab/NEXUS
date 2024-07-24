#include <seal/seal.h>

#include <chrono>
#include <iostream>
#include <vector>

#include "Bootstrapper.h"
#include "ckks_evaluator.h"

using namespace std;
using namespace chrono;
using namespace seal;

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

  // CKKS
  KeyGenerator *keygen = nullptr;

 public:
  CKKSEvaluator *ckks = nullptr;
  GaloisKeys *bootstrapping_keys = nullptr;

  // execution time minus bootstrapping setup time
  chrono::duration<int64_t, nano> time_elapsed = chrono::duration<int64_t, nano>(0);

  ArgmaxEvaluator(CKKSEvaluator &ckks, KeyGenerator &keygen, int main_mod_count) {
    this->ckks = &ckks;
    this->keygen = &keygen;

    GaloisKeys boot_keys;
    this->bootstrapping_keys = &boot_keys;

    this->logN = log2(ckks.degree);
    this->logn = logN - 2;
    this->sparse_slots = (1 << logn);

    this->total_level = main_mod_count + bs_mod_count;
  }

  void argmax(Ciphertext &x, Ciphertext &res, int len);

  void bootstrap(Ciphertext &x);
};
