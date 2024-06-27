#include "softmax.cuh"

using namespace nexus;

void SoftmaxEvaluator::softmax2(PhantomCiphertext &x, PhantomCiphertext &res, int len) {
  PhantomCiphertext tmp, a, b, sign, a_plus_b, a_minus_b, exp_x;
  PhantomPlaintext zero_point_five, one, inverse_64;
  int log_step = log2(len);

  ckks->encoder.encode(ckks->init_vec_with_value(ckks->slot_count, 1.0 / 64), x.params_id(), x.scale(), inverse_64);
  ckks->evaluator.multiply_plain_inplace(x, inverse_64);
  ckks->evaluator.rescale_to_next_inplace(x);
  ckks->encoder.encode(ckks->init_vec_with_value(ckks->slot_count, 1.0), x.params_id(), x.scale(), one);
  ckks->evaluator.add_plain_inplace(x, one);

  // x^64
  ckks->evaluator.square(x, x);
  ckks->evaluator.relinearize_inplace(x, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(x);
  ckks->evaluator.square(x, x);
  ckks->evaluator.relinearize_inplace(x, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(x);
  ckks->evaluator.square(x, x);
  ckks->evaluator.relinearize_inplace(x, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(x);
  ckks->evaluator.square(x, x);
  ckks->evaluator.relinearize_inplace(x, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(x);
  ckks->evaluator.square(x, x);
  ckks->evaluator.relinearize_inplace(x, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(x);
  ckks->evaluator.square(x, x);
  ckks->evaluator.relinearize_inplace(x, *ckks->relin_keys);
  ckks->evaluator.rescale_to_next_inplace(x);

  res = x;
}
