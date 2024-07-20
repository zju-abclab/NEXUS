#include "ModularReducer.h"

ModularReducer::ModularReducer(long _boundary_K, double _log_width, long _deg, long _num_double_formula, long _inverse_deg,
                               SEALContext &_context,
                               CKKSEncoder &_encoder,
                               Encryptor &_encryptor,
                               Evaluator &_evaluator,
                               RelinKeys &_relin_keys,
                               Decryptor &_decryptor) : boundary_K(_boundary_K), log_width(_log_width), deg(_deg), num_double_formula(_num_double_formula), inverse_deg(_inverse_deg), context(_context), encoder(_encoder), encryptor(_encryptor), evaluator(_evaluator), relin_keys(_relin_keys), decryptor(_decryptor) {
  inverse_log_width = -log2(sin(2 * M_PI * pow(2.0, -log_width)));
  poly_generator = new RemezCos(rmparm, boundary_K, log_width, deg, (1 << num_double_formula));
  inverse_poly_generator = new RemezArcsin(rmparm, inverse_log_width, inverse_deg);

  inverse_poly_generator->params.log_scan_step_diff = 12;
  inverse_poly_generator->params.RR_prec = 1000;
}

void ModularReducer::double_angle_formula(Ciphertext &cipher) {
  evaluator.square_inplace(cipher);
  evaluator.relinearize_inplace(cipher, relin_keys);
  evaluator.rescale_to_next_inplace(cipher);
  evaluator.double_inplace(cipher);
  evaluator.add_const(cipher, -1.0, cipher);
}

void ModularReducer::double_angle_formula_scaled(Ciphertext &cipher, double scale_coeff) {
  evaluator.square_inplace(cipher);
  evaluator.relinearize_inplace(cipher, relin_keys);
  evaluator.rescale_to_next_inplace(cipher);
  evaluator.double_inplace(cipher);
  evaluator.add_const(cipher, -scale_coeff, cipher);
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

void ModularReducer::modular_reduction(Ciphertext &rtn, Ciphertext &cipher) {
  Ciphertext tmp1, tmp2;
  Plaintext tmpplain;
  tmp1 = cipher;

  sin_cos_polynomial.homomorphic_poly_evaluation(context, encoder, encryptor, evaluator, relin_keys, tmp2, tmp1, decryptor);

  if (inverse_deg == 1) {
    double curr_scale = scale_inverse_coeff;
    for (int i = 0; i < num_double_formula; i++) {
      curr_scale = curr_scale * curr_scale;
      double_angle_formula_scaled(tmp2, curr_scale);
    }
    rtn = tmp2;
  }

  else {
    for (int i = 0; i < num_double_formula; i++) double_angle_formula(tmp2);
    inverse_sin_polynomial.homomorphic_poly_evaluation(context, encoder, encryptor, evaluator, relin_keys, rtn, tmp2, decryptor);
  }
}
