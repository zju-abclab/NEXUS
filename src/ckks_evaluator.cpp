#include "ckks_evaluator.h"

#include <chrono>
#include <iostream>

using namespace chrono;

// @deprecated Bootstrapping has been implemented
void CKKSEvaluator::re_encrypt(Ciphertext &ct) {
  auto start = high_resolution_clock::now();
  while (ct.coeff_modulus_size() > 1) {
    evaluator->mod_switch_to_next_inplace(ct);
  }

  vector<seal_byte> data;
  data.resize(ct.save_size(compr_mode_type::zstd));
  comm += ct.save(data.data(), data.size(), compr_mode_type::zstd);
  // cout << "Communication cost:  " <<
  // ct.save(data.data(),data.size(),compr_mode_type::zstd) << " bytes" << endl;

  Plaintext temp;
  vector<double> v;
  decryptor->decrypt(ct, temp);
  encoder->decode(temp, v);

  encoder->encode(v, ct.scale(), temp);
  encryptor->encrypt(temp, ct);

  auto end = high_resolution_clock::now();
  cout << duration_cast<milliseconds>(end - start).count() / 2 << " milliseconds" << endl;

  // cout << "depth = " <<
  // context->get_context_data(ct.parms_id())->chain_index() << "\n";
}

void CKKSEvaluator::print_decrypted_ct(Ciphertext &ct, int nums) {
  Plaintext temp;
  vector<double> v;
  decryptor->decrypt(ct, temp);
  encoder->decode(temp, v);
  for (int i = 0; i < nums; i++) {
    cout << v[i] << " ";
  }
  cout << "\n";
}

void CKKSEvaluator::print_decoded_pt(Plaintext &pt, int num) {
  vector<double> v;

  encoder->decode(pt, v);

  for (int i = 0; i < num; i++) {
    cout << v[i] << " ";
  }
  cout << endl;
}

double CKKSEvaluator::calculateMAE(vector<double> &y_true, Ciphertext &ct, int N) {
  Plaintext temp;
  vector<double> y_pred;
  decryptor->decrypt(ct, temp);
  encoder->decode(temp, y_pred);

  double sum_absolute_errors = 0.0;
  for (size_t i = 0; i < N; ++i) {
    sum_absolute_errors += abs(y_true[i] - y_pred[i]);
  }

  return sum_absolute_errors / N;
}

vector<double> CKKSEvaluator::init_vec_with_value(int N, double init_value) {
  std::vector<double> v(N);

  for (int i = 0; i < N; ++i) {
    v[i] = init_value;
  }

  return v;
}

vector<double> CKKSEvaluator::init_mask(int N, int m) {
  std::vector<double> v(N);

  for (int i = 0; i < N / m; ++i) {
    if (i % 2 == 0) {
      for (int j = 0; j < m; ++j) {
        v[i * m + j] = 1;
      }
    } else {
      for (int j = 0; j < m; ++j) {
        v[i * m + j] = 0;
      }
    }
  }

  return v;
}

Ciphertext CKKSEvaluator::poly_eval(Ciphertext x, vector<Plaintext> coeff) {
  // cout << "Initial depth " <<
  // context->get_context_data(x.parms_id())->chain_index() << "\n"; x^2

  Ciphertext x_2;
  evaluator->square(x, x_2);
  evaluator->relinearize_inplace(x_2, *relin_keys);
  evaluator->rescale_to_next_inplace(x_2);
  // cout << "x: " << x.scale() << " x_2: " << x_2.scale() << "\n";

  // x^3
  Ciphertext x_3;
  evaluator->mod_switch_to_inplace(x, x_2.parms_id());
  evaluator->multiply(x_2, x, x_3);
  evaluator->relinearize_inplace(x_3, *relin_keys);
  evaluator->rescale_to_next_inplace(x_3);
  // cout << "x_2: " << x_2.scale() << " x_3: " << x_3.scale() << "\n";

  // x^4
  Ciphertext x_4;
  evaluator->square(x_2, x_4);
  evaluator->relinearize_inplace(x_4, *relin_keys);
  evaluator->rescale_to_next_inplace(x_4);

  // x^5
  Ciphertext x_5;

  // x^7
  Ciphertext x_7;

  Ciphertext new_x, new_x3, sum_5_7;

  evaluator->mod_switch_to_inplace(coeff[1], x.parms_id());
  evaluator->multiply_plain(x, coeff[1], new_x);
  evaluator->rescale_to_next_inplace(new_x);

  evaluator->mod_switch_to_inplace(coeff[3], x_3.parms_id());
  evaluator->multiply_plain(x_3, coeff[3], new_x3);
  evaluator->rescale_to_next_inplace(new_x3);

  // cout << x.scale() << " " << x_3.scale() << " " << scale << endl;

  evaluator->mod_switch_to_inplace(coeff[5], x_3.parms_id());
  evaluator->mod_switch_to_inplace(x, x_3.parms_id());
  evaluator->multiply_plain(x, coeff[5], x_5);
  evaluator->rescale_to_next_inplace(x_5);

  evaluator->mod_switch_to_inplace(coeff[7], x_3.parms_id());
  evaluator->multiply_plain(x_3, coeff[7], x_7);
  evaluator->rescale_to_next_inplace(x_7);

  // cout << x_5.scale() << " " << x_7.scale() << " " << scale << endl;

  x_5.scale() = scale;
  x_7.scale() = scale;

  evaluator->add(x_5, x_7, sum_5_7);

  evaluator->mod_switch_to_next_inplace(x_4);
  evaluator->multiply_inplace(sum_5_7, x_4);
  evaluator->relinearize_inplace(sum_5_7, *relin_keys);
  evaluator->rescale_to_next_inplace(sum_5_7);

  // c7*x^7 + c5*x^5 + c3*x^3 + c1*x
  new_x.scale() = scale;
  new_x3.scale() = scale;
  sum_5_7.scale() = scale;

  evaluator->mod_switch_to_inplace(new_x, sum_5_7.parms_id());
  evaluator->mod_switch_to_inplace(new_x3, sum_5_7.parms_id());

  evaluator->add_inplace(sum_5_7, new_x);
  evaluator->add_inplace(sum_5_7, new_x3);

  // cout << "Final depth " <<
  // context->get_context_data(x.parms_id())->chain_index() << "\n";

  return sum_5_7;
}

// Output range: [-sgn_factor, sgn_factor]
Ciphertext CKKSEvaluator::sgn_eval(Ciphertext x, int d_g, int d_f, double sgn_factor) {
  // Compute sign function coefficients
  vector<double> f4_coeffs = F4_COEFFS;
  vector<double> g4_coeffs = G4_COEFFS;
  vector<double> f4_coeffs_last(10, 0.0);
  vector<double> g4_coeffs_last(10, 0.0);

  for (int i = 0; i <= 9; i++) {
    f4_coeffs[i] /= f4_scale;
    f4_coeffs_last[i] = f4_coeffs[i] * sgn_factor;

    g4_coeffs[i] /= g4_scale;
    g4_coeffs_last[i] = g4_coeffs[i] * sgn_factor;
  }

  Ciphertext dest = x;
  for (int i = 0; i < d_g; i++) {
    if (i == d_g - 1) {
      eval_odd_deg9_poly(g4_coeffs_last, dest, dest);
    } else {
      eval_odd_deg9_poly(g4_coeffs, dest, dest);
    }
  }
  for (int i = 0; i < d_f; i++) {
    if (i == d_f - 1) {
      eval_odd_deg9_poly(f4_coeffs_last, dest, dest);
    } else {
      eval_odd_deg9_poly(f4_coeffs, dest, dest);
    }
  }
  return dest;
}

void CKKSEvaluator::eval_odd_deg9_poly(vector<double> &a, Ciphertext &x, Ciphertext &dest) {
  /*
      (polyeval/odd9.h)
      P(x) = a9 x^9 + a7 x^7 + a5 x^5 + a3 x^3 + a1 x

      T1 = (a3 + a5 x^2) x^3
      T2 = (a7 x + a9 x^3) x^6
      T3 = a1 x
      P(x) = T1 + T2 + T3

      Depth=4, #Muls=5

      Exactly what babystep_giantstep would do, but written explicitly to optimize

      ###

      -> Errorless Polynomial Evaluation (3.2. of https://eprint.iacr.org/2020/1203)
      GOAL: evaluate a polynomial exactly so no need to stabilize and lose precision
      (x at level L and scale D --> P(x) at level L-4 and scale D)
      it's possible to do this exactly for polyeval as (x,x2,x3,x6) determine the scale D_L for each involved level L:
      (assume the primes at levels L to L-4 are p, q, r, s)

      level       ctx       scale (D_l)
      ==================================
        L          x          D
        L-1        x2         D^2 / p
        L-2        x3         D^3 / pq
        L-3        x6         D^6 / p^2 q^2 r

      Just by encoding constants at different scales we can make every ctx at level l be at scale D_l
      (not possible in general, e.g. rescale(x2*x2) produces L-2 ciphertext with scale D^4/ppq)
      (to fix this we would use the Adjust op. that multiplies ctx by constants and Algo 3 for primes from
     https://eprint.iacr.org/2020/1118)

      Now we know that sc(P(x)) should be D, so we recursively go back to compute the scales for each coefficient
      sc(T1)=sc(T2)=sc(T3)=sc(P(x))=D

      T3:
          sc(a1) = q (should be p but it gets multiplied with modswitched x)

      T2:
          sc(x^6) = D^6 / p^2 q^2 r, so sc(a7*x) = sc(a9*x^3) = p^2 q^2 r s / D^5
          next, sc(a7) = p^2 q^3 r s / D^6
          similarly, sc(a9) = p^3 q^3 r^2 s / D^8

      T1:
          sc(x^3) = D^3 / pq
          implying sc(a3) = pqr / D^2 and also sc(a5*x^2) = pqr / D^2
          as sc(x^2) = D^2 / p this implies sc(a5) = p^2 q^2 r / D^4
  */
  // chrono::high_resolution_clock::time_point time_start, time_end;
  // time_start = high_resolution_clock::now();
  double D = x.scale();  // maybe not init_scale but preserved

  uint64_t p = get_modulus(x, 1);
  uint64_t q = get_modulus(x, 2);
  uint64_t r = get_modulus(x, 3);
  uint64_t s = get_modulus(x, 4);
  uint64_t t = get_modulus(x, 5);

  p = q;
  q = r;
  r = s;
  s = t;

  vector<double> a_scales(10);
  a_scales[1] = q;
  a_scales[3] = (double)p / D * q / D * r;
  a_scales[5] = (double)p / D * p / D * q / D * q / D * r;
  a_scales[7] = (double)p / D * p / D * q / D * q / D * q / D * r / D * s;
  a_scales[9] = (double)p / D * p / D * p / D * q / D * q / D * q / D * r / D * r / D * s;

  ///////////////////////////////////////////////
  Ciphertext x2, x3, x6;

  evaluator->square(x, x2);
  evaluator->relinearize_inplace(x2, *relin_keys);
  evaluator->rescale_to_next_inplace(x2);  // L-1

  evaluator->mod_switch_to_next_inplace(x);  // L-1
  evaluator->multiply(x2, x, x3);
  evaluator->relinearize_inplace(x3, *relin_keys);
  evaluator->rescale_to_next_inplace(x3);  // L-2

  evaluator->square(x3, x6);
  evaluator->relinearize_inplace(x6, *relin_keys);
  evaluator->rescale_to_next_inplace(x6);  // L-3

  Plaintext a1, a3, a5, a7, a9;

  // Build T1
  Ciphertext T1;
  double a5_scale = D / x2.scale() * p / x3.scale() * q;
  encoder->encode(a[5], x2.parms_id(), a5_scale, a5);  // L-1
  evaluator->multiply_plain(x2, a5, T1);
  evaluator->rescale_to_next_inplace(T1);  // L-2

  // Update: using a_scales[3] is only approx. correct, so we directly use T1.scale()
  encoder->encode(a[3], T1.parms_id(), T1.scale(), a3);  // L-2

  evaluator->add_plain_inplace(T1, a3);  // L-2

  evaluator->multiply_inplace(T1, x3);
  evaluator->relinearize_inplace(T1, *relin_keys);
  evaluator->rescale_to_next_inplace(T1);  // L-3

  // Build T2
  Ciphertext T2;
  Plaintext a9_switched;
  double a9_scale = D / x3.scale() * r / x6.scale() * q;
  encoder->encode(a[9], x3.parms_id(), a9_scale, a9);  // L-2
  evaluator->multiply_plain(x3, a9, T2);
  evaluator->rescale_to_next_inplace(T2);  // L-3

  Ciphertext a7x;
  double a7_scale = T2.scale() / x.scale() * p;
  encoder->encode(a[7], x.parms_id(), a7_scale, a7);  // L-1 (x was modswitched)
  evaluator->multiply_plain(x, a7, a7x);
  evaluator->rescale_to_next_inplace(a7x);               // L-2
  evaluator->mod_switch_to_inplace(a7x, T2.parms_id());  // L-3

  double mid_scale = (T2.scale() + a7x.scale()) / 2;
  T2.scale() = a7x.scale() = mid_scale;  // this is the correct scale now, need to set it still to avoid SEAL assert
  evaluator->add_inplace(T2, a7x);       // L-3
  evaluator->multiply_inplace(T2, x6);
  evaluator->relinearize_inplace(T2, *relin_keys);
  evaluator->rescale_to_next_inplace(T2);  // L-4

  // Build T3
  Ciphertext T3;
  encoder->encode(a[1], x.parms_id(), p, a1);  // L-1 (x was modswitched)
  evaluator->multiply_plain(x, a1, T3);
  evaluator->rescale_to_next_inplace(T3);  // L-2

  // T1, T2 and T3 should be on the same scale up to floating point
  // but we still need to set them manually to avoid SEAL assert
  double mid3_scale = (T1.scale() + T2.scale() + T3.scale()) / 3;
  T1.scale() = T2.scale() = T3.scale() = mid3_scale;

  dest = T2;
  evaluator->mod_switch_to_inplace(T1, dest.parms_id());  // L-4
  evaluator->add_inplace(dest, T1);
  evaluator->mod_switch_to_inplace(T3, dest.parms_id());  // L-4
  evaluator->add_inplace(dest, T3);

  /////////////////////////////////////////
  // it should be ==D but we don't stabilize if it's not, D' != D is ok
  // the goal was to make T1+T2+T3 work with minimal loss in precision
  // time_end = high_resolution_clock::now();
  // cout << "Poly eval took " << duration_cast<milliseconds>(time_end - time_start).count() << " ms" << endl;
}

Ciphertext CKKSEvaluator::newtonIter(Ciphertext x, Ciphertext res, int iter) {
  for (int i = 0; i < iter; i++) {
    // if (context->get_context_data(res.parms_id())->chain_index() < 4)
    //     re_encrypt(res);
    // cout << i << " " << depth(res) << "\n";
    Plaintext three_half, neg_half;
    encoder->encode(1.5, scale, three_half);
    encoder->encode(-0.5, scale, neg_half);

    // x^2
    Ciphertext res_sq;
    evaluator->square(res, res_sq);
    evaluator->relinearize_inplace(res_sq, *relin_keys);
    evaluator->rescale_to_next_inplace(res_sq);
    // printVector(res_sq);
    // evaluator->negate_inplace(res_sq);
    // printVector(res_sq, 3);
    // cout << "square\n";

    //-0.5*x*b
    Ciphertext res_x;
    evaluator->mod_switch_to_inplace(neg_half, x.parms_id());
    evaluator->multiply_plain(x, neg_half, res_x);
    evaluator->rescale_to_next_inplace(res_x);
    if (context->get_context_data(res.parms_id())->chain_index() <
        context->get_context_data(res_x.parms_id())->chain_index())
      evaluator->mod_switch_to_inplace(res_x, res.parms_id());
    else
      evaluator->mod_switch_to_inplace(res, res_x.parms_id());

    evaluator->multiply_inplace(res_x, res);
    evaluator->relinearize_inplace(res_x, *relin_keys);
    evaluator->rescale_to_next_inplace(res_x);
    // cout << "negate\n";

    //-0.5*b*x^3
    evaluator->mod_switch_to_inplace(res_sq, res_x.parms_id());
    evaluator->multiply_inplace(res_x, res_sq);
    evaluator->relinearize_inplace(res_x, *relin_keys);
    evaluator->rescale_to_next_inplace(res_x);
    // cout << "res_x\n";
    // printVector(res_x, 3);
    // 1.5*x
    evaluator->mod_switch_to_inplace(three_half, res.parms_id());
    evaluator->multiply_plain_inplace(res, three_half);
    evaluator->rescale_to_next_inplace(res);
    // cout << "constant\n";

    //-0.5*b*x^3 + 1.5*x
    evaluator->mod_switch_to_inplace(res, res_x.parms_id());
    res_x.scale() = scale;
    res.scale() = scale;
    evaluator->add_inplace(res, res_x);
    // cout << "final\n";
  }
  return res;
}

pair<Ciphertext, Ciphertext> CKKSEvaluator::goldSchmidtIter(Ciphertext v, Ciphertext y, int d) {
  Ciphertext x, h, r, temp;
  Plaintext constant;
  encoder->encode(0.5, scale, constant);

  // GoldSchmidt's algorithm
  evaluator->mod_switch_to_inplace(v, y.parms_id());
  evaluator->multiply(v, y, x);
  evaluator->relinearize_inplace(x, *relin_keys);
  evaluator->rescale_to_next_inplace(x);
  evaluator->mod_switch_to_inplace(constant, y.parms_id());
  evaluator->multiply_plain(y, constant, h);
  evaluator->rescale_to_next_inplace(h);

  for (int i = 0; i < d; i++) {
    encoder->encode(0.5, scale, constant);

    evaluator->multiply(x, h, r);
    evaluator->relinearize_inplace(r, *relin_keys);
    evaluator->rescale_to_next_inplace(r);
    r.scale() = scale;
    evaluator->negate(r, temp);
    evaluator->mod_switch_to_inplace(constant, temp.parms_id());
    evaluator->add_plain(temp, constant, r);
    // cout << "r\n";

    // x = x + x*r
    evaluator->mod_switch_to_inplace(x, r.parms_id());
    evaluator->multiply(x, r, temp);
    evaluator->relinearize_inplace(temp, *relin_keys);
    evaluator->rescale_to_next_inplace(temp);
    x.scale() = scale;
    temp.scale() = scale;
    evaluator->mod_switch_to_inplace(x, temp.parms_id());
    evaluator->add_inplace(x, temp);
    // cout << "x\n";

    // h = h + h*r
    evaluator->mod_switch_to_inplace(h, r.parms_id());
    evaluator->multiply(h, r, temp);
    evaluator->relinearize_inplace(temp, *relin_keys);
    evaluator->rescale_to_next_inplace(temp);
    h.scale() = scale;
    temp.scale() = scale;
    evaluator->mod_switch_to_inplace(h, temp.parms_id());
    evaluator->add_inplace(h, temp);
    // cout << "h\n";
  }
  encoder->encode(2.0, scale, constant);
  evaluator->mod_switch_to_inplace(constant, h.parms_id());
  evaluator->multiply_plain_inplace(h, constant);
  evaluator->rescale_to_next_inplace(h);

  return make_pair(x, h);
}

Ciphertext CKKSEvaluator::invert_sqrt(Ciphertext x, int d_newt, int d_gold) {
  Ciphertext res = initGuess(x);
  Ciphertext y = newtonIter(x, res, d_newt);
  pair<Ciphertext, Ciphertext> sqrt_inv_sqrt = goldSchmidtIter(x, y, d_gold);
  return sqrt_inv_sqrt.second;
}

Ciphertext CKKSEvaluator::inverse(Ciphertext x, int iter) {
  Ciphertext y, tmp, res;
  Plaintext one;
  encoder->encode(1.0, x.parms_id(), x.scale(), one);
  evaluator->sub_plain(x, one, y);
  evaluator->negate_inplace(y);
  evaluator->add_plain(y, one, tmp);
  res = tmp;
  for (int i = 0; i < iter; i++) {
    evaluator->square_inplace(y);
    evaluator->relinearize_inplace(y, *relin_keys);
    evaluator->rescale_to_next_inplace(y);

    encoder->encode(1.0, y.parms_id(), y.scale(), one);
    evaluator->add_plain(y, one, tmp);

    evaluator->mod_switch_to_inplace(res, tmp.parms_id());
    evaluator->multiply_inplace(res, tmp);
    evaluator->relinearize_inplace(res, *relin_keys);
    evaluator->rescale_to_next_inplace(res);
  }
  return res;
}

uint64_t CKKSEvaluator::get_modulus(Ciphertext &x, int k) {
  const vector<Modulus> &modulus = context->get_context_data(x.parms_id())->parms().coeff_modulus();
  int sz = modulus.size();
  return modulus[sz - k].value();
}

Ciphertext CKKSEvaluator::initGuess(Ciphertext x) {
  Plaintext A, B;
  encoder->encode(-1.29054537e-04, scale, A);
  encoder->encode(1.29054537e-01, scale, B);
  return evalLine(x, A, B);
}

Ciphertext CKKSEvaluator::evalLine(Ciphertext x, Plaintext m, Plaintext c) {
  evaluator->mod_switch_to_inplace(m, x.parms_id());
  evaluator->multiply_plain_inplace(x, m);
  evaluator->rescale_to_next_inplace(x);
  evaluator->mod_switch_to_inplace(c, x.parms_id());
  x.scale() = scale;
  evaluator->add_plain_inplace(x, c);
  return x;
}

Ciphertext CKKSEvaluator::exp(Ciphertext x) {
  Plaintext one, inverse_128;
  encoder->encode(0.0078125, x.parms_id(), x.scale(), inverse_128);
  evaluator->multiply_plain_inplace(x, inverse_128);
  evaluator->rescale_to_next_inplace(x);
  encoder->encode(1.0, x.parms_id(), x.scale(), one);
  evaluator->add_plain_inplace(x, one);
  // x^128
  for (int i = 0; i < log2(128); i++) {
    evaluator->square(x, x);
    evaluator->relinearize_inplace(x, *relin_keys);
    evaluator->rescale_to_next_inplace(x);
  }
  return x;
}
