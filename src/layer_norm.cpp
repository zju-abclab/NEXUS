#include "layer_norm.h"

#include <iostream>
#include <vector>

void LNEvaluator::layer_norm(Ciphertext &x, Ciphertext &res, int len) {
  Ciphertext tmp, x2;
  int log_step = log2(len);
  ckks->evaluator->rotate_vector(x, -len, *ckks->galois_keys, tmp);
  ckks->evaluator->add_inplace(x, tmp);
  ckks->evaluator->square(x, x2);
  ckks->evaluator->relinearize_inplace(x2, *ckks->relin_keys);
  ckks->evaluator->rescale_to_next_inplace(x2);
  tmp = x2;

  for (int i = 0; i < log_step; ++i) {
    ckks->evaluator->rotate_vector(tmp, pow(2, i), *ckks->galois_keys, res);
    ckks->evaluator->add_inplace(res, tmp);
    tmp = res;
  }

  Plaintext delta;
  ckks->encoder->encode(1.0 / 768, res.parms_id(), res.scale(), delta);
  ckks->evaluator->multiply_plain_inplace(res, delta);
  ckks->evaluator->rescale_to_next_inplace(res);

  res = ckks->invert_sqrt(res, 4, 2);

  ckks->evaluator->mod_switch_to_inplace(x, res.parms_id());
  ckks->evaluator->multiply(res, x, res);
  ckks->evaluator->relinearize_inplace(res, *ckks->relin_keys);
  ckks->evaluator->rescale_to_next_inplace(res);
}
