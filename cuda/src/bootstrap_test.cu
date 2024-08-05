#include <chrono>
#include <random>

#include "bootstrapping/Bootstrapper.cuh"
#include "phantom.h"

using namespace std;
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

int main() {
  long boundary_K = 25;
  long deg = 59;
  long scale_factor = 2;
  long inverse_deg = 1;

  // The following parameters have been adjusted to satisfy the memory constraints of an A100 GPU
  long logN = 15;  // 16 -> 15
  long loge = 10;

  long logn = 13;  // 14 -> 13
  long sparse_slots = (1 << logn);

  int logp = 46;
  int logq = 51;
  int log_special_prime = 51;

  int secret_key_hamming_weight = 192;

  int remaining_level = 16;
  int boot_level = 14;  // Greater than: subsum 1 + coefftoslot 2 + ModReduction 9 + slottocoeff 2
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

  std::cout << "Setting Parameters..." << endl;
  phantom::EncryptionParameters parms(scheme_type::ckks);
  size_t poly_modulus_degree = (size_t)(1 << logN);
  double scale = pow(2.0, logp);

  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(phantom::arith::CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec));
  parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
  parms.set_sparse_slots(sparse_slots);

  PhantomContext context(parms);

  PhantomSecretKey secret_key(context);
  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);
  PhantomGaloisKey galois_keys;

  PhantomCKKSEncoder encoder(context);

  CKKSEvaluator ckks_evaluator(&context, &public_key, &secret_key, &encoder, &relin_keys, &galois_keys, scale);

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

  std::cout << "Generating Optimal Minimax Polynomials..." << endl;
  bootstrapper.prepare_mod_polynomial();

  std::cout << "Adding Bootstrapping Keys..." << endl;
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logN - 1; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  bootstrapper.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);

  ckks_evaluator.decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks_evaluator.galois_keys));
  std::cout << "Galois key generated from steps vector." << endl;

  bootstrapper.slot_vec.push_back(logn);

  std::cout << "Generating Linear Transformation Coefficients..." << endl;
  bootstrapper.generate_LT_coefficient_3();

  vector<double> sparse(sparse_slots, 0.0);
  vector<double> input(slot_count, 0.0);
  vector<double> before(slot_count, 0.0);
  vector<double> after(slot_count, 0.0);

  random_real(sparse, sparse_slots);

  PhantomPlaintext plain;
  PhantomCiphertext cipher;

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
  bootstrapper.bootstrap_3(rtn, cipher);

  duration<double> sec = system_clock::now() - start;
  std::cout << "Bootstrapping took: " << sec.count() << "s" << endl;
  std::cout << "Return cipher level: " << rtn.coeff_modulus_size() << endl;

  ckks_evaluator.decryptor.decrypt(rtn, plain);
  ckks_evaluator.encoder.decode(plain, after);

  double mean_err = 0;
  for (long i = 0; i < sparse_slots; i++) {
    if (i < 10) std::cout << before[i] << " <----> " << after[i] << endl;
    mean_err += abs(before[i] - after[i]);
  }
  mean_err /= sparse_slots;
  std::cout << "Mean absolute error: " << mean_err << endl;
}
