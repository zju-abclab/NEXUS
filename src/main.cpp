#include <seal/seal.h>
#include <seal/util/uintarith.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "argmax.h"
#include "gelu.h"
#include "layer_norm.h"
#include "matrix_mul.h"
#include "softmax.h"

using namespace std;
using namespace seal;
using namespace seal::util;
using namespace std::chrono;

void MM_test();
void argmax_test();

// Choose a test target here:
int TEST_TARGET_IDX = 0;

vector<string> TEST_TARGETS = {"MatMul", "Argmax", "GELU", "LayerNorm", "SoftMax"};
vector<int> COEFF_MODULI = {58, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 58};
vector<int> MM_COEFF_MODULI = {60, 40, 60};
double SCALE = pow(2.0, 40);

string TEST_TARGET = TEST_TARGETS[TEST_TARGET_IDX];

int main() {
  if (TEST_TARGET == TEST_TARGETS[0]) {
    MM_test();
    return 0;
  }

  if (TEST_TARGET == TEST_TARGETS[1]) {
    argmax_test();
    return 0;
  }

  long logN = 16;
  size_t poly_modulus_degree = 1 << logN;

  EncryptionParameters parms(scheme_type::ckks);
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, COEFF_MODULI));

  SEALContext context(parms, true, sec_level_type::none);

  KeyGenerator keygen(context);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  keygen.create_public_key(public_key);
  RelinKeys relin_keys;
  keygen.create_relin_keys(relin_keys);
  GaloisKeys galois_keys;
  keygen.create_galois_keys(galois_keys);

  Encryptor encryptor(context, public_key);
  CKKSEncoder encoder(context);
  Evaluator evaluator(context, encoder);
  Decryptor decryptor(context, secret_key);

  CKKSEvaluator ckks_evaluator(context, encryptor, decryptor, encoder, evaluator, SCALE, relin_keys, galois_keys);
  GeLUEvaluator gelu_evaluator(ckks_evaluator);
  LNEvaluator ln_evaluator(ckks_evaluator);
  SoftmaxEvaluator softmax_evaluator(ckks_evaluator);

  Plaintext plain_input;
  Ciphertext cipher_input;
  Ciphertext cipher_output;
  vector<double> output;

  /*
      GELU
  */
  if (TEST_TARGET == TEST_TARGETS[2]) {
    double num;
    vector<double> input, gelu_calibration;

    ifstream input_file("../data/input/gelu_input_32768.txt");
    while (input_file >> num) {
      input.push_back(num);
    }
    input_file.close();

    ifstream calibration_file("../data/calibration/gelu_calibration_32768.txt");
    while (calibration_file >> num) {
      gelu_calibration.push_back(num);
    }
    calibration_file.close();

    ckks_evaluator.encoder->encode(input, SCALE, plain_input);
    ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);

    auto start = high_resolution_clock::now();
    gelu_evaluator.gelu(cipher_input, cipher_output);
    auto end = high_resolution_clock::now();

    cout << "[GELU] 32768 takes: " << duration_cast<milliseconds>(end - start).count() << " milliseconds" << endl;
    cout << "Mean Absolute Error: "
         << ckks_evaluator.calculateMAE(gelu_calibration, cipher_output, poly_modulus_degree / 2)
         << endl;

    return 0;
  }

  /*
      LayerNorm
  */
  if (TEST_TARGET == TEST_TARGETS[3]) {
    double num;
    vector<double> input, layernorm_calibration;

    ifstream input_file("../data/input/layernorm_input_16_768.txt");
    while (input_file >> num) {
      input.push_back(num);
    }
    input_file.close();

    ifstream calibration_file("../data/calibration/layernorm_calibration_16_768.txt");
    while (calibration_file >> num) {
      layernorm_calibration.push_back(num);
    }
    calibration_file.close();

    ckks_evaluator.encoder->encode(input, SCALE, plain_input);
    ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);

    auto start = high_resolution_clock::now();
    ln_evaluator.layer_norm(cipher_input, cipher_output, 1024);
    auto end = high_resolution_clock::now();

    cout << "[LayerNorm] 16 x 768 takes: " << duration_cast<milliseconds>(end - start).count() << " milliseconds" << endl;
    cout << "Mean Absolute Error: "
         << ckks_evaluator.calculateMAE(layernorm_calibration, cipher_output, 768)
         << endl;

    return 0;
  }

  /*
      Softmax
  */
  if (TEST_TARGET == TEST_TARGETS[4]) {
    double num;
    vector<double> input, softmax_calibration;

    ifstream input_file("../data/input/softmax_input_128_128.txt");
    while (input_file >> num) {
      input.push_back(num);
    }
    input_file.close();

    ifstream calibration_file("../data/calibration/softmax_calibration_128_128.txt");
    while (calibration_file >> num) {
      softmax_calibration.push_back(num);
    }
    calibration_file.close();

    ckks_evaluator.encoder->encode(input, SCALE, plain_input);
    ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);

    auto start = high_resolution_clock::now();
    softmax_evaluator.softmax(cipher_input, cipher_output, 128);
    auto end = high_resolution_clock::now();

    cout << "[Softmax] 128 x 128 takes: " << duration_cast<milliseconds>(end - start).count() << " milliseconds"
         << endl;
    cout << "Mean Absolute Error: " << ckks_evaluator.calculateMAE(softmax_calibration, cipher_output, 128)
         << endl;

    return 0;
  }
}

void argmax_test() {
  long logN = 15;
  long logn = logN - 2;
  long sparse_slots = (1 << logn);

  int logp = 46;
  int logq = 51;
  int log_special_prime = 58;

  // QuickMax: 17
  int main_mod_count = 17;  // Mod count after bootstrapping: 18
  // Subsum 1 + coefftoslot 2 + ModReduction 9 + slottocoeff 2
  int bs_mod_count = 14;
  int total_level = main_mod_count + bs_mod_count;

  // Bootstrapping parameters
  long boundary_K = 25;
  long deg = 59;
  long scale_factor = 2;
  long inverse_deg = 1;
  long loge = 10;

  int secret_key_hamming_weight = 192;

  vector<int> coeff_bit_vec;
  coeff_bit_vec.push_back(logq);
  for (int i = 0; i < main_mod_count; i++) {
    coeff_bit_vec.push_back(logp);
  }
  for (int i = 0; i < bs_mod_count; i++) {
    coeff_bit_vec.push_back(logq);
  }
  coeff_bit_vec.push_back(log_special_prime);

  size_t poly_modulus_degree = (size_t)(1 << logN);
  double scale = pow(2.0, logp);

  EncryptionParameters parms(scheme_type::ckks);
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec));
  parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
  parms.set_sparse_slots(sparse_slots);

  SEALContext context(parms, true, sec_level_type::none);

  KeyGenerator keygen(context);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  keygen.create_public_key(public_key);
  RelinKeys relin_keys;
  keygen.create_relin_keys(relin_keys);
  GaloisKeys galois_keys;
  keygen.create_galois_keys(galois_keys);
  GaloisKeys bootstrapping_keys;

  CKKSEncoder encoder(context);
  Encryptor encryptor(context, public_key);
  Evaluator evaluator(context, encoder);
  Decryptor decryptor(context, secret_key);

  CKKSEvaluator ckks_evaluator(context, encryptor, decryptor, encoder, evaluator, scale, relin_keys, galois_keys);
  Bootstrapper bootstrapper(loge, logn, logN - 1, total_level, scale, boundary_K, deg, scale_factor, inverse_deg,
                            context, keygen, encoder, encryptor, decryptor, evaluator,
                            relin_keys, bootstrapping_keys);

  cout << "Generating Optimal Minimax Polynomials..." << endl;
  bootstrapper.prepare_mod_polynomial();

  cout << "Adding Bootstrapping Keys..." << endl;
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logN - 1; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  bootstrapper.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  keygen.create_galois_keys(gal_steps_vector, bootstrapping_keys);
  bootstrapper.slot_vec.push_back(logn);

  cout << "Generating Linear Transformation Coefficients..." << endl;
  bootstrapper.generate_LT_coefficient_3();

  ArgmaxEvaluator argmax_evaluator(ckks_evaluator, bootstrapper);

  size_t slot_count = encoder.slot_count();

  Plaintext plain_input;
  Ciphertext cipher_input;
  Ciphertext cipher_output;
  vector<double> output;
  vector<double> input(slot_count, 0.0);

  double num;
  int argmax_input_size = 0;
  vector<double> argmax_input(sparse_slots, 0.0), argmax_calibration;

  ifstream input_file("../data/input/argmax_input_8.txt");
  while (input_file >> num) {
    argmax_input[argmax_input_size] = num;
    argmax_input_size++;
  }
  input_file.close();

  ifstream calibration_file("../data/calibration/argmax_calibration_8.txt");
  while (calibration_file >> num) {
    argmax_calibration.push_back(num);
  }
  calibration_file.close();

  // Spare input (TODO: create a dedicated encoding function: encode_sparse in ckks evaluator)
  for (size_t i = 0; i < slot_count; i++) {
    input[i] = argmax_input[i % sparse_slots];
  }

  ckks_evaluator.encoder->encode(input, scale, plain_input);
  ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);

  // Mod switch to remaining level
  for (int i = 0; i < bs_mod_count; i++) {
    ckks_evaluator.evaluator->mod_switch_to_next_inplace(cipher_input);
  }

  auto start = high_resolution_clock::now();
  argmax_evaluator.argmax(cipher_input, cipher_output, argmax_input_size);
  auto end = high_resolution_clock::now();

  cout << "[Argmax] 32768 takes: " << duration_cast<milliseconds>(end - start).count() << " milliseconds"
       << endl;
  cout << "Mean Absolute Error: " << ckks_evaluator.calculateMAE(argmax_calibration, cipher_output, argmax_input_size) << endl;
}

void MM_test() {
  long logN = 13;
  size_t poly_modulus_degree = 1 << logN;

  EncryptionParameters parms(scheme_type::ckks);
  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, MM_COEFF_MODULI));

  SEALContext context(parms, true, sec_level_type::none);

  KeyGenerator keygen(context);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  keygen.create_public_key(public_key);
  RelinKeys relin_keys;
  keygen.create_relin_keys(relin_keys);
  GaloisKeys galois_keys;

  std::vector<std::uint32_t> rots;
  for (int i = 0; i < logN; i++) {
    rots.push_back((poly_modulus_degree + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
  }
  keygen.create_galois_keys(rots, galois_keys);

  Encryptor encryptor(context, public_key);
  CKKSEncoder encoder(context);
  Evaluator evaluator(context, encoder);
  Decryptor decryptor(context, secret_key);

  CKKSEvaluator ckks_evaluator(context, encryptor, decryptor, encoder, evaluator, SCALE, relin_keys, galois_keys);
  MMEvaluator mme(ckks_evaluator);

  std::vector<std::vector<double>> matrix_4096x768 = mme.readMatrix("../data/input/matrixmul_input_m_128_n_768_k_64_batch_128.txt", 4096, 768);
  std::vector<std::vector<double>> matrix_768x64 = mme.readMatrix("../data/input/matrix_input_n_768_k_64.txt", 768, 64);

  auto matrix_4096x768_T = mme.transposeMatrix(matrix_4096x768);
  auto matrix_768x64_T = mme.transposeMatrix(matrix_768x64);

  std::vector<std::vector<double>> row_pack;
  std::vector<double> row_ct(poly_modulus_degree, 0.0);
  for (auto i = 0; i < 64 * 768; i++) {
    int row = i / 768;
    int col = i % 768;
    row_ct[i % poly_modulus_degree] = matrix_768x64_T[row][col];
    if (i % poly_modulus_degree == (poly_modulus_degree - 1)) {
      row_pack.push_back(row_ct);
    }
  }

  vector<Ciphertext> res;
  auto start = high_resolution_clock::now();
  mme.matrix_mul(matrix_4096x768_T, row_pack, res);
  auto end = high_resolution_clock::now();
  cout << "[MatMul] 4096x768 x 768x64 takes: " << duration_cast<milliseconds>(end - start).count() << " milliseconds"
       << endl;

  std::vector<std::vector<double>> matrix_4096x64 = mme.readMatrix("../data/calibration/matrix_output_m_128_k_64_batch_128.txt", 4096, 64);
  auto matrix_4096x64_T = mme.transposeMatrix(matrix_4096x64);

  // Calculate the error of the first column
  double average_err = 0.0;
  Plaintext res_pt;
  vector<double> mm_res;
  ckks_evaluator.decryptor->decrypt(res[0], res_pt);
  ckks_evaluator.encoder->decode(res_pt, mm_res);

  for (auto i = 0; i < 4096; i++) {
    average_err += fabs(mm_res[i] / 2.0 - matrix_4096x64_T[0][i]);
  }

  std::cout << "Average error: " << average_err / 4096.0 << std::endl;
}
