#include "ModularReducer.cuh"

ModularReducer::ModularReducer(long _boundary_K, double _log_width, long _deg, long _num_double_formula, long _inverse_deg,
                               CKKSEvaluator *_ckks) : boundary_K(_boundary_K), log_width(_log_width), deg(_deg), num_double_formula(_num_double_formula), inverse_deg(_inverse_deg), ckks(_ckks) {
  inverse_log_width = -log2(sin(2 * M_PI * pow(2.0, -log_width)));
  poly_generator = new RemezCos(rmparm, boundary_K, log_width, deg, (1 << num_double_formula));
  inverse_poly_generator = new RemezArcsin(rmparm, inverse_log_width, inverse_deg);

  inverse_poly_generator->params.log_scan_step_diff = 12;
  inverse_poly_generator->params.RR_prec = 1000;
}

void ModularReducer::double_angle_formula(PhantomCiphertext &cipher) {
  ckks->evaluator.square_inplace(cipher);
  ckks->evaluator.relinearize_inplace(cipher, *(ckks->relin_keys));
  ckks->evaluator.rescale_to_next_inplace(cipher);
  ckks->evaluator.double_inplace(cipher);
  ckks->evaluator.add_const(cipher, -1.0, cipher);
}

void ModularReducer::double_angle_formula_scaled(PhantomCiphertext &cipher, double scale_coeff) {
  ckks->evaluator.square_inplace(cipher);
  ckks->evaluator.relinearize_inplace(cipher, *(ckks->relin_keys));
  ckks->evaluator.rescale_to_next_inplace(cipher);
  ckks->evaluator.double_inplace(cipher);
  ckks->evaluator.add_const(cipher, -scale_coeff, cipher);
}

void ModularReducer::generate_sin_cos_polynomial() {
  poly_generator->generate_optimal_poly(sin_cos_polynomial);
  sin_cos_polynomial.generate_poly_heap();
}

void ModularReducer::generate_inverse_sine_polynomial() {
  inverse_poly_generator->generate_optimal_poly(inverse_sin_polynomial);
  if (inverse_deg > 3) inverse_sin_polynomial.generate_poly_heap_odd();
  if (inverse_deg == 1) {
    scale_inverse_coeff = to_double(inverse_sin_polynomial.coeff[1]);
    for (int i = 0; i < num_double_formula; i++) scale_inverse_coeff = sqrt(scale_inverse_coeff);
    sin_cos_polynomial.constmul(to_RR(scale_inverse_coeff));
    sin_cos_polynomial.generate_poly_heap();
  }
}

void ModularReducer::write_polynomials() {
  ofstream sin_cos_out("cosine.txt"), inverse_out("inverse_sine.txt");
  sin_cos_polynomial.write_heap_to_file(sin_cos_out);
  inverse_sin_polynomial.write_heap_to_file(inverse_out);
  sin_cos_out.close();
  inverse_out.close();
}

void ModularReducer::modular_reduction(PhantomCiphertext &rtn, PhantomCiphertext &cipher) {
  PhantomCiphertext tmp1, tmp2;
  PhantomPlaintext tmpplain;
  tmp1 = cipher;

  sin_cos_polynomial.homomorphic_poly_evaluation(ckks, tmp2, tmp1);

  // ckks->print_decrypted_ct(tmp2, 10);

  if (inverse_deg == 1) {
    double curr_scale = scale_inverse_coeff;
    for (int i = 0; i < num_double_formula; i++) {
      curr_scale = curr_scale * curr_scale;
      double_angle_formula_scaled(tmp2, curr_scale);
    }
    rtn = tmp2;
  } else {
    for (int i = 0; i < num_double_formula; i++) double_angle_formula(tmp2);
    inverse_sin_polynomial.homomorphic_poly_evaluation(ckks, rtn, tmp2);
  }
}
