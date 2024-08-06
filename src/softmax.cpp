#include "softmax.h"

#include <iostream>
#include <vector>

void SoftmaxEvaluator::softmax(Ciphertext &x, Ciphertext &res, int len) {
  Ciphertext tmp, exp_x;

  int log_step = log2(len);

  ckks->evaluator->rotate_vector(x, -len, *ckks->galois_keys, tmp);
  ckks->evaluator->add_inplace(x, tmp);

  exp_x = ckks->exp(x);

  tmp = exp_x;
  for (int i = 0; i < log_step; ++i) {
    ckks->evaluator->rotate_vector(tmp, pow(2, i), *ckks->galois_keys, res);
    ckks->evaluator->add_inplace(res, tmp);
    tmp = res;
  }

  // let res/delta in [0, 1]
  Plaintext delta;
  ckks->encoder->encode(0.01, res.parms_id(), res.scale(), delta);
  ckks->evaluator->multiply_plain_inplace(res, delta);
  ckks->evaluator->rescale_to_next_inplace(res);

  res = ckks->inverse(res);

  // recover to 1/res
  ckks->encoder->encode(0.01, res.parms_id(), res.scale(), delta);
  ckks->evaluator->multiply_plain_inplace(res, delta);
  ckks->evaluator->rescale_to_next_inplace(res);

  ckks->evaluator->mod_switch_to_inplace(exp_x, res.parms_id());
  ckks->evaluator->multiply(res, exp_x, res);
  ckks->evaluator->relinearize_inplace(res, *ckks->relin_keys);
  ckks->evaluator->rescale_to_next_inplace(res);

  // cout << "Moduli left after SoftMax: " << res.coeff_modulus_size() << endl;
}
