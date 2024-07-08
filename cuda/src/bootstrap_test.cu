#include <NTL/RR.h>

#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>

#include "bootstrapping/Bootstrapper.cuh"
#include "bootstrapping/ModularReducer.cuh"
#include "bootstrapping/common/Polynomial.cuh"
#include "gelu.cuh"
#include "layer_norm.cuh"
#include "matrix_mul.cuh"
#include "phantom.h"
#include "softmax.cuh"

using namespace std;
using namespace NTL;
using namespace phantom;
using namespace chrono;

void random_real(vector<double> &vec, size_t size) {
  random_device rn;
  mt19937_64 rnd(rn());
  thread_local std::uniform_real_distribution<double> distribution(-1, 1);

  vec.reserve(size);

  for (size_t i = 0; i < size; i++) {
    vec[i] = distribution(rnd);
  }
}

double recover_vefore(double before, long boundary_K) {
  double scaled = boundary_K * before;
  return scaled - round(scaled);
}

int main() {
  long boundary_K = 25;
  long deg = 59;
  long scale_factor = 2;
  long inverse_deg = 1;

  long logN = 15;
  long loge = 10;

  long logn = 14;
  long logn_2 = 13;
  long logn_3 = 12;
  long sparse_slots = (1 << logn);

  int logp = 46;
  int logq = 51;
  int log_special_prime = 51;

  int secret_key_hamming_weight = 192;

  // int remaining_level = 14; // Calculation required
  int remaining_level = 16;  // Calculation required
  int boot_level = 14;       // greater than: subsum 1 + coefftoslot 2 + ModReduction 9 + slottocoeff 2
  int total_level = remaining_level + boot_level;

  vector<int> coeff_bit_vec;
  coeff_bit_vec.push_back(logq);
  for (int i = 0; i < remaining_level; i++) {
    coeff_bit_vec.push_back(logp);
  }

  for (int i = 0; i < boot_level; i++) {
    coeff_bit_vec.push_back(logq);
  }
  coeff_bit_vec.push_back(log_special_prime);

  cout << "Setting Parameters" << endl;
  EncryptionParameters parms(scheme_type::ckks);
  size_t poly_modulus_degree = (size_t)(1 << logN);
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec));
  double scale = pow(2.0, logp);
  // modified SEAL
  parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
  // parms.set_sparse_slots(sparse_slots);

  PhantomContext context(parms);

  PhantomSecretKey secret_key(context);
  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);
  PhantomGaloisKey galois_keys;

  PhantomCKKSEncoder encoder(context);

  CKKSEvaluator ckks_evaluator(context, public_key, secret_key, encoder, relin_keys, galois_keys, scale);

  size_t slot_count = encoder.slot_count();

  Bootstrapper bootstrapper(
      loge,
      logn,
      logN - 1,
      total_level,
      scale,
      boundary_K,
      deg,
      scale_factor,
      inverse_deg,
      ckks_evaluator);
  // Bootstrapper bootstrapper_2(
  //     loge,
  //     logn_2,
  //     logN - 1,
  //     total_level,
  //     scale,
  //     boundary_K,
  //     deg,
  //     scale_factor,
  //     inverse_deg,
  //     ckks_evaluator);
  // Bootstrapper bootstrapper_3(
  //     loge,
  //     logn_3,
  //     logN - 1,
  //     total_level,
  //     scale,
  //     boundary_K,
  //     deg,
  //     scale_factor,
  //     inverse_deg,
  //     ckks_evaluator);

  cout << "Generating Optimal Minimax Polynomials..." << endl;
  cout << "-1" << endl;
  bootstrapper.prepare_mod_polynomial();
  // bootstrapper_2.prepare_mod_polynomial();
  // bootstrapper_3.prepare_mod_polynomial();
  cout << "Adding Bootstrapping Keys..." << endl;
  // bootstrapper.addBootKeys_3_other_slots(gal_keys, slot_vec);
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logN - 1; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  bootstrapper.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  cout << "-2" << endl;
  // bootstrapper_2.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  // bootstrapper_3.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  ckks_evaluator.decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks_evaluator.galois_keys));
  cout << "-3" << endl;

  bootstrapper.slot_vec.push_back(logn);
  // bootstrapper_2.slot_vec.push_back(logn_2);
  // bootstrapper_3.slot_vec.push_back(logn_3);

  cout << "Generating Linear Transformation Coefficients..." << endl;
  bootstrapper.generate_LT_coefficient_3();
  cout << "-4" << endl;
  // bootstrapper_2.generate_LT_coefficient_3();
  // bootstrapper_3.generate_LT_coefficient_3();

  double tot_err = 0, mean_err;
  size_t iterations = 1;

  vector<double> sparse(sparse_slots);
  vector<double> input(slot_count);
  vector<double> before(slot_count);
  vector<double> after(slot_count);

  random_real(sparse, sparse_slots);

  PhantomPlaintext plain;
  PhantomCiphertext cipher;

  for (size_t _ = 0; _ < iterations; _++) {
    if (_ == 0)
      sparse_slots = (1 << logn);
    else if (_ == 1)
      sparse_slots = (1 << logn_2);
    else if (_ == 2)
      sparse_slots = (1 << logn_3);

    cout << _ << "-th iteration : sparse_slots = " << sparse_slots << endl;

    for (size_t i = 0; i < slot_count; i++) {
      input[i] = sparse[i % sparse_slots];
    }

    ckks_evaluator.encoder.encode(input, scale, plain);
    ckks_evaluator.encryptor.encrypt(plain, cipher);

    for (int i = 0; i < total_level; i++) {
      ckks_evaluator.evaluator.mod_switch_to_next_inplace(cipher);
    }

    PhantomCiphertext rtn;

    ckks_evaluator.decryptor.decrypt(cipher, plain);
    // encoder.decode(plain, before, sparse_slots);
    ckks_evaluator.encoder.decode(plain, before);

    auto start = system_clock::now();

    if (_ == 0)
      bootstrapper.bootstrap_3(rtn, cipher);
    // else if (_ == 1)
      // bootstrapper_2.bootstrap_3(rtn, cipher);
    // else if (_ == 2)
      // bootstrapper_3.bootstrap_3(rtn, cipher);

    duration<double> sec = system_clock::now() - start;
    cout << "bootstrapping time : " << sec.count() << "s" << endl;

    cout << "return cipher level: " << rtn.coeff_modulus_size() << endl;

    rtn = ckks_evaluator.sgn_eval(rtn, 2, 2);

    ckks_evaluator.decryptor.decrypt(rtn, plain);
    // encoder.decode(plain, after, sparse_slots);
    ckks_evaluator.encoder.decode(plain, after);

    for (long i = 0; i < sparse_slots; i++) {
      if (before[i] > 0) {
        before[i] = 0.5;
      } else {
        before[i] = -0.5;
      }
    }

    mean_err = 0;
    for (long i = 0; i < sparse_slots; i++) {
      // cout << i << "m: " << recover_vefore(before[i].real(), boundary_K) << "\td: " << after[i].real()<< endl;

      if (i < 10)
        cout << i << " " << before[i] << "<---->" << after[i] << endl;

      mean_err += abs(before[i] - after[i]);
      // if (file.is_open())
      // {
      //     file << before[i].real() - after[i].real() << ","
      //      << before[i].imag() - after[i].imag() << "," << flush;
      // }
    }
    mean_err /= sparse_slots;
    cout << "Absolute mean of error: " << mean_err << endl;
    tot_err += mean_err;
  }
  tot_err /= iterations;

  // if (file.is_open())
  // {
  //     file << endl;
  // }
  cout << " mean error: " << tot_err << endl;
  // file.close();

  return 0;
}
