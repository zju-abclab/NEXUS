#include "layer_norm.cuh"

using namespace nexus;

void LNEvaluator::layer_norm(PhantomCiphertext &a, PhantomCiphertext &y, int len) {
  PhantomCiphertext tmp, x2;

  int log_step = log2(len);
  ckks->evaluator.rotate_vector(a, -len, *ckks->galois_keys, tmp);
  ckks->evaluator.add_inplace(a, tmp);

  ckks->evaluator.square(a, x2);
  ckks->evaluator.relinearize_inplace(x2, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(x2);

  tmp = x2;
  for (int i = 0; i < log_step; ++i) {
    ckks->evaluator.rotate_vector(tmp, pow(2, i), *ckks->galois_keys, y);
    ckks->evaluator.add_inplace(y, tmp);
    tmp = y;
  }

  PhantomPlaintext delta;
  ckks->encoder.encode(1.0 / 768, y.params_id(), y.scale(), delta);
  ckks->evaluator.multiply_plain_inplace(y, delta);
  ckks->evaluator.rescale_to_next_inplace(y);

  y = ckks->invert_sqrt(y, 4, 2);

  ckks->evaluator.mod_switch_to_inplace(a, y.params_id());
  ckks->evaluator.multiply(y, a, y);
  ckks->evaluator.relinearize_inplace(y, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(y);

  // cout << "Moduli left after LayerNorm: " << y.coeff_modulus_size() << endl;
}
