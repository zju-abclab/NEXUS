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

  long logN = 15;  // 16 -> 15
  long loge = 10;

  long logn = 13;  // 14 -> 13
  // long logn_2 = 13;
  // long logn_3 = 12;
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

  cout << "Setting Parameters..." << endl;
  phantom::EncryptionParameters parms(scheme_type::ckks);
  size_t poly_modulus_degree = (size_t)(1 << logN);
  double scale = pow(2.0, logp);

  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(phantom::arith::CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec));
  parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
  parms.set_sparse_slots(sparse_slots);

  PhantomContext context(parms);

  PhantomSecretKey secret_key(context);
  std::ifstream sk_in("../../sk.txt");
  secret_key.load_secret_key(context, sk_in);

  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);
  PhantomGaloisKey galois_keys;

  // Bootstrapping evaluation function test Galois Key
  // PhantomGaloisKey galois_keys = secret_key.create_galois_keys(context);

  PhantomCKKSEncoder encoder(context);

  CKKSEvaluator ckks_evaluator(&context, &public_key, &secret_key, &encoder, &relin_keys, &galois_keys, scale);

  // Bootstrapping evaluation function tests ------------------------------------------
  // vector<double> output;
  // vector<complex<double>> output2;
  // vector<double> input = {1.0, 2.0, 3.0, 4.0, 5.0};
  // vector<double> input2 = {1.0, 2.0, 3.0, 4.0};
  // vector<complex<double>> input3 = {1.0 + 2.0i, 3.0 + 4.0i, 5.0 + 6.0i};
  // PhantomPlaintext plain;
  // PhantomPlaintext plain2;
  // PhantomPlaintext plain3;
  // PhantomPlaintext output_plain;
  // PhantomPlaintext output2_plain;
  // PhantomCiphertext cipher;
  // PhantomCiphertext cipher2;
  // PhantomCiphertext cipher3;
  // PhantomCiphertext output_cipher;
  // ckks_evaluator.encoder.encode(input, scale, plain);
  // // ckks_evaluator.encoder.encode(input2, scale, plain2);
  // // ckks_evaluator.encoder.encode(input3, scale, plain3);
  // ckks_evaluator.encryptor.encrypt(plain, cipher);
  // // ckks_evaluator.encryptor.encrypt(plain2, cipher2);
  // // ckks_evaluator.encryptor.encrypt(plain3, cipher3);

  // auto cipher_data = new uint64_t[cipher.poly_modulus_degree() * cipher.coeff_modulus_size()];
  // cudaMemcpy(cipher_data, cipher.data(), cipher.poly_modulus_degree() * cipher.coeff_modulus_size() * sizeof(uint64_t), cudaMemcpyDeviceToHost);
  // for (int i = 0; i < 10; i++) {
  //   cout << cipher_data[i] << " ";
  // }
  // cout << endl;

  // ckks_evaluator.evaluator.add_inplace_reduced_error(cipher, cipher2);
  // ckks_evaluator.evaluator.sub_inplace_reduced_error(cipher, cipher2);

  // ckks_evaluator.evaluator.multiply_inplace_reduced_error(cipher, cipher2, *ckks_evaluator.relin_keys);
  // ckks_evaluator.evaluator.multiply_vector_inplace_reduced_error(cipher, input2);

  // ckks_evaluator.evaluator.add_const_inplace(cipher, 1.0);
  // ckks_evaluator.evaluator.multiply_const_inplace(cipher, 2.0);
  // ckks_evaluator.evaluator.add_const_inplace(cipher, -1.0);

  // ckks_evaluator.evaluator.complex_conjugate_inplace(cipher, *ckks_evaluator.galois_keys);
  // ckks_evaluator.evaluator.rotate_vector_inplace(cipher, 1, *ckks_evaluator.galois_keys);

  // ckks_evaluator.evaluator.multiply_vector_inplace_reduced_error(cipher, input3);

  // cout << "Cipher:" << endl;
  // cout << "chain_index: " << context.get_context_data(cipher.params_id()).chain_depth() << endl;
  // cout << "coeff_modulus_size: " << cipher.coeff_modulus_size() << endl;
  // cout << "scale: " << cipher.scale() << endl;

  // ckks_evaluator.decryptor.decrypt(cipher, output_plain);
  // ckks_evaluator.encoder.decode(output_plain, output);

  // for (int i = 0; i < 5; i++) {
  //   cout << output[i] << " ";
  // }
  // cout << endl;
  // Bootstrapping evaluation function tests ------------------------------------------

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
      &ckks_evaluator);
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
  bootstrapper.prepare_mod_polynomial();  // Good to go
  // bootstrapper_2.prepare_mod_polynomial();
  // bootstrapper_3.prepare_mod_polynomial();

  cout << "Generating Galois Keys..." << endl;
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logN - 1; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  bootstrapper.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  // bootstrapper_2.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  // bootstrapper_3.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);

  ckks_evaluator.decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks_evaluator.galois_keys));  // Good to go
  cout << "Galois key generated from steps vector." << endl;

  bootstrapper.slot_vec.push_back(logn);
  // bootstrapper_2.slot_vec.push_back(logn_2);
  // bootstrapper_3.slot_vec.push_back(logn_3);

  cout << "Generating Linear Transformation Coefficients..." << endl;
  bootstrapper.generate_LT_coefficient_3();  // Good to go
  // bootstrapper_2.generate_LT_coefficient_3();
  // bootstrapper_3.generate_LT_coefficient_3();

  double tot_err = 0, mean_err;
  size_t iterations = 1;

  vector<double> sparse(sparse_slots, 0.0);
  vector<double> input(slot_count, 0.0);
  vector<double> before(slot_count, 0.0);
  vector<double> after(slot_count, 0.0);

  random_real(sparse, sparse_slots);

  PhantomPlaintext plain;
  PhantomCiphertext cipher;

  for (size_t _ = 0; _ < iterations; _++) {
    if (_ == 0)
      sparse_slots = (1 << logn);
    // else if (_ == 1)
    // sparse_slots = (1 << logn_2);
    // else if (_ == 2)
    // sparse_slots = (1 << logn_3);

    cout << _ << "-th iteration : sparse_slots = " << sparse_slots << endl;

    // Create input cipher
    for (size_t i = 0; i < slot_count; i++) {
      input[i] = sparse[i % sparse_slots];
    }

    ckks_evaluator.encoder.encode(input, scale, plain);
    ckks_evaluator.encryptor.encrypt(plain, cipher);

    // Mod switch to the lowest level
    for (int i = 0; i < total_level; i++) {
      ckks_evaluator.evaluator.mod_switch_to_next_inplace(cipher);
    }

    // Decrypt input cipher to obtain the original input
    ckks_evaluator.decryptor.decrypt(cipher, plain);
    ckks_evaluator.encoder.decode(plain, before);

    auto start = system_clock::now();

    PhantomCiphertext rtn;

    if (_ == 0)
      bootstrapper.bootstrap_3(rtn, cipher);
    // else if (_ == 1)
    // bootstrapper_2.bootstrap_3(rtn, cipher);
    // else if (_ == 2)
    // bootstrapper_3.bootstrap_3(rtn, cipher);

    duration<double> sec = system_clock::now() - start;
    cout << "Bootstrapping took: " << sec.count() << "s" << endl;
    cout << "Return cipher level: " << rtn.coeff_modulus_size() << endl;

    // rtn = ckks_evaluator.sgn_eval(rtn, 2, 2);

    ckks_evaluator.decryptor.decrypt(rtn, plain);
    ckks_evaluator.encoder.decode(plain, after);

    // for (long i = 0; i < sparse_slots; i++) {
    //   if (before[i] > 0) {
    //     before[i] = 0.5;
    //   } else {
    //     before[i] = -0.5;
    //   }
    // }

    mean_err = 0;
    for (long i = 0; i < sparse_slots; i++) {
      if (i < 10)
        cout << before[i] << " <----> " << after[i] << endl;
      mean_err += abs(before[i] - after[i]);
    }
    mean_err /= sparse_slots;
    cout << "Mean absolute error: " << mean_err << endl;
    tot_err += mean_err;
  }

  tot_err /= iterations;
  cout << "Mean error: " << tot_err << endl;

  return 0;
}
