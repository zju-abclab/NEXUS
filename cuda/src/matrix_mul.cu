#include "matrix_mul.cuh"

using namespace nexus;

vector<PhantomCiphertext> MMEvaluator::expand_ciphertext(
    const PhantomCiphertext &encrypted, uint32_t m, PhantomGaloisKey &galkey, vector<uint32_t> &galois_elts) {
  uint32_t logm = ceil(log2(m));
  auto n = ckks->degree;

  vector<PhantomCiphertext> temp;
  temp.push_back(encrypted);

  PhantomCiphertext tempctxt;
  PhantomCiphertext tempctxt_rotated;
  PhantomCiphertext tempctxt_shifted;
  PhantomCiphertext tempctxt_rotatedshifted;

  for (uint32_t i = 0; i < logm; i++) {
    vector<PhantomCiphertext> newtemp(temp.size() << 1);
    int index_raw = (n << 1) - (1 << i);
    int index = (index_raw * galois_elts[i]) % (n << 1);
    for (uint32_t a = 0; a < temp.size(); a++) {
      ckks->evaluator.apply_galois(temp[a], ckks->rots[i], *(ckks->galois_keys), tempctxt_rotated);  // sub
      ckks->evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
      multiply_power_of_x(temp[a], tempctxt_shifted, index_raw);  // x**-1
      multiply_power_of_x(tempctxt_rotated, tempctxt_rotatedshifted, index);
      ckks->evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
    }
    temp = newtemp;
  }
  return temp;
}
