#include "argmax.h"

#include <chrono>
#include <iostream>
#include <vector>

using namespace chrono;

void ArgmaxEvaluator::argmax(Ciphertext &x, Ciphertext &res, int len) {
  // TODO: Assert len is << N, and a power of 2

  Ciphertext tmp, a, b, sign, a_plus_b, a_minus_b, a_minus_b_sgn;
  Plaintext one, half;

  int log_step = log2(len);

  // Transform x = [a_0, ..., a_n, 0, ..., 0] to [a_0, ..., a_n, a_0, ..., a_n, 0, ..., 0]
  ckks->evaluator->rotate_vector(x, -len, *ckks->galois_keys, tmp);
  ckks->evaluator->add_inplace(x, tmp);

  a = x;
  for (int i = 0; i < log_step; ++i) {
    ckks->evaluator->rotate_vector(a, pow(2, i), *ckks->galois_keys, b);

    ckks->evaluator->add(a, b, a_plus_b);
    ckks->evaluator->sub(a, b, a_minus_b);
    sign = ckks->sgn_eval(a_minus_b, 2, 2);

    // (a - b) * sgn(a - b) / 2
    ckks->evaluator->mod_switch_to_inplace(a_minus_b, sign.parms_id());
    ckks->evaluator->multiply(a_minus_b, sign, a_minus_b_sgn);
    ckks->evaluator->relinearize_inplace(a_minus_b_sgn, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(a_minus_b_sgn);

    // (a + b) / 2
    ckks->encoder->encode(0.5, a_plus_b.parms_id(), a_plus_b.scale(), half);
    ckks->evaluator->multiply_plain_inplace(a_plus_b, half);
    ckks->evaluator->rescale_to_next_inplace(a_plus_b);

    // a = max(a, b)
    a_plus_b.scale() = a_minus_b_sgn.scale();
    ckks->evaluator->mod_switch_to_inplace(a_plus_b, a_minus_b_sgn.parms_id());
    ckks->evaluator->add(a_plus_b, a_minus_b_sgn, a);

    bootstrap(a);
  }

  x.scale() = a.scale();
  ckks->evaluator->mod_switch_to_inplace(x, a.parms_id());
  ckks->evaluator->sub(x, a, res);

  res = ckks->sgn_eval(res, 2, 2, 1.0);

  ckks->encoder->encode(1.0, res.parms_id(), res.scale(), one);
  ckks->evaluator->add_plain_inplace(res, one);
}

void ArgmaxEvaluator::bootstrap(Ciphertext &x) {
  cout << "Bootstrapping started" << endl;

  if (x.coeff_modulus_size() > 1) {
    cout << "Ciphertext is not at lowest level, remaining level(s): " + to_string(x.coeff_modulus_size()) << endl;
    cout << "Mod switching to the lowest level..." << endl;

    while (x.coeff_modulus_size() > 1) {
      ckks->evaluator->mod_switch_to_next_inplace(x);
    }
  }

  Ciphertext x_0 = x;
  bootstrapper->set_final_scale(x.scale());
  cout << "Bootstrapping..." << endl;

  auto start = system_clock::now();

  bootstrapper->bootstrap_3(x, x_0);

  duration<double> sec = system_clock::now() - start;
  cout << "Bootstrapping took: " << sec.count() << "s" << endl;
  cout << "New ciphertext depth: " << x.coeff_modulus_size() << endl;
}
