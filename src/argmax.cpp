#include "argmax.h"

void ArgmaxEvaluator::argmax(Ciphertext &x, Ciphertext &res, int len) {
  // TODO: Assert len is << N, and a power of 2

  auto start = high_resolution_clock::now();

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

    auto end = high_resolution_clock::now();
    time_elapsed += (end - start);
    // ckks->print_decrypted_ct(a, 8);
    bootstrap(a);
    // ckks->print_decrypted_ct(a, 8);
    start = high_resolution_clock::now();
  }

  // ckks->print_decrypted_ct(a, 8);

  x.scale() = a.scale();
  ckks->evaluator->mod_switch_to_inplace(x, a.parms_id());
  ckks->evaluator->sub(x, a, res);

  res = ckks->sgn_eval(res, 2, 2, 1.0);

  // ckks->print_decrypted_ct(a, 8);

  ckks->encoder->encode(1.0, res.parms_id(), res.scale(), one);
  ckks->evaluator->add_plain_inplace(res, one);

  auto end = high_resolution_clock::now();
  time_elapsed += (end - start);
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
  Bootstrapper bootstrapper(
      loge,
      logn,
      logN - 1,
      total_level,
      ckks->scale,
      boundary_K,
      deg,
      scale_factor,
      inverse_deg,
      *ckks->context,
      *keygen,
      *ckks->encoder,
      *ckks->encryptor,
      *ckks->decryptor,
      *ckks->evaluator,
      *ckks->relin_keys,
      *bootstrapping_keys);

  cout << "Generating Optimal Minimax Polynomials..." << endl;
  bootstrapper.prepare_mod_polynomial();

  cout << "Adding Bootstrapping Keys..." << endl;
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logN - 1; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  bootstrapper.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  keygen->create_galois_keys(gal_steps_vector, *bootstrapping_keys);
  bootstrapper.slot_vec.push_back(logn);

  cout << "Generating Linear Transformation Coefficients..." << endl;
  bootstrapper.generate_LT_coefficient_3();

  cout << "Bootstrapping..." << endl;
  auto start = system_clock::now();
  bootstrapper.bootstrap_3(x, x_0);
  duration<double> sec = system_clock::now() - start;
  cout << "Bootstrapping took: " << sec.count() << "s" << endl;
  cout << "New ciphertext depth: " << x.coeff_modulus_size() << endl;
}
